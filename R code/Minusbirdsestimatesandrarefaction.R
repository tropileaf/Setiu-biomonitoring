library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(iNEXT)
library(ggplot2)
library(writexl)
library(lubridate)

#read and clean detection data
animals <- read_excel("Combined_nohuman.xlsx",
                      sheet = "-birds"
) %>%
  mutate(
    Method = str_trim(as.character(Method)),
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT-S" ~ "CT_S",
      TRUE ~ Method
    ),
    Site = str_trim(as.character(Site)),
    Camera_ID = str_trim(as.character(Camera_ID)),
    
    # fixes A1cam -> A1, S1cam -> S1 if needed
    Camera_ID = str_replace(Camera_ID, "cam$", ""),
    
    Date = as.Date(`Day/Month`),
    Visit_No = as.integer(`Visit No.`),
    PC_No = as.integer(`PC No.`),
    Taxon = str_squish(paste(Genus, Species))
  ) %>%
  filter(
    Method %in% c("CT_A", "CT_S", "PC"),
    !is.na(Taxon),
    Taxon != "NA NA",
    Taxon != ""
  )
animals_clean <- animals %>%
  mutate(
    Taxon = trimws(as.character(Taxon)),
    Taxon = na_if(Taxon, ""),
    Location = trimws(as.character(Location)),
    Location = na_if(Location, ""),
    Site = trimws(as.character(Site)),
    Site = na_if(Site, ""),
    Camera_ID = trimws(as.character(Camera_ID)),
    Camera_ID = na_if(Camera_ID, ""),
    Visit_No = as.character(Visit_No),
    PC_No = as.character(PC_No),
    Date = as.character(Date)
  )

#create PC sampling units 
pc_detected <- animals %>%
  filter(Method == "PC") %>%
  mutate(
    sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")
  ) %>%
  distinct(sample_id, Site, Visit_No, PC_No, Taxon)

pc_full <- expand_grid(
  Site = unique(animals$Site[animals$Method == "PC"]),
  Visit_No = 1:3,
  PC_No = 1:10
) %>%
  mutate(
    Method = "PC",
    sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")
  )

pc_all <- pc_full %>%
  left_join(pc_detected, by = c("sample_id", "Site", "Visit_No", "PC_No"))

#double check effort PC
pc_all %>%
  distinct(Site, sample_id) %>%
  count(Site)

#create CT days
deploy <- read_excel("camera_deployment_dates.xlsx") %>%
  mutate(
    Method = str_trim(as.character(Method)),
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT-S" ~ "CT_S",
      TRUE ~ Method
    ),
    Site = str_trim(as.character(Site)),
    Camera_ID = str_trim(as.character(Camera_ID)),
    Camera_ID = str_replace(Camera_ID, "cam$", ""),
    start_date = as.Date(start_date),
    end_date = as.Date(end_date)
  )
#generate active days
ct_full <- deploy %>%
  rowwise() %>%
  mutate(Date = list(seq(start_date, end_date, by = "day"))) %>%
  unnest(Date) %>%
  ungroup() %>%
  mutate(
    sample_id = paste(Method, Site, Camera_ID, Date, sep = "_")
  )
#CT detections
ct_detected <- animals %>%
  filter(Method %in% c("CT_A", "CT_S")) %>%
  mutate(
    sample_id = paste(Method, Site, Camera_ID, Date, sep = "_")
  ) %>%
  distinct(Method, Site, Camera_ID, Date, sample_id, Taxon)
#join detections and active days
ct_all <- ct_full %>%
  left_join(
    ct_detected,
    by = c("Method", "Site", "Camera_ID", "Date", "sample_id")
  )
#check effort
ct_all %>%
  distinct(Method, Site, Camera_ID, sample_id) %>%
  count(Method, Site, Camera_ID)

ct_all %>%
  count(Method, is.na(Taxon))

#build PC iNEXT input
pc_incidence <- pc_all %>%
  filter(!is.na(Taxon)) %>%
  distinct(sample_id, Taxon) %>%
  group_by(Taxon) %>%
  summarise(incidence = n(), .groups = "drop")

pc_inext <- c(
  nrow(pc_full),
  pc_incidence$incidence
)

#build CT iNEXT
ct_effort <- ct_all %>%
  distinct(Method, sample_id) %>%
  count(Method, name = "n_sampling_units")

ct_incidence <- ct_all %>%
  filter(!is.na(Taxon)) %>%
  distinct(Method, sample_id, Taxon) %>%
  group_by(Method, Taxon) %>%
  summarise(incidence = n(), .groups = "drop")

ct_input <- ct_incidence %>%
  left_join(ct_effort, by = "Method") %>%
  group_by(Method) %>%
  summarise(
    input = list(c(first(n_sampling_units), incidence)),
    .groups = "drop"
  )

ct_list <- setNames(ct_input$input, ct_input$Method)

#combine
bird_data <- c(
  ct_list,
  list(PC = pc_inext)
)

names(bird_data)
sapply(bird_data, function(x) x[1])

#run standardized to 30 days
out_30 <- iNEXT(
  bird_data,
  q = 0,
  datatype = "incidence_freq",
  endpoint = 30
)
#plot 
p_30 <- ggiNEXT(out_30, type = 1) +
  theme_classic() +
  labs(
    x = "Number of sampling units",
    y = "Species richness",
    color = "Method",
    fill = "Method"
  ) +
  theme(
    legend.title = element_blank()
  ) +
  guides(shape = "none")

p_30

#full effort
out_full <- iNEXT(
  bird_data,
  q = 0,
  datatype = "incidence_freq"
)

p_full <- ggiNEXT(out_full, type = 1) +
  theme_classic(base_size = 18) +
  labs(
    x = "Number of sampling units",
    y = "Species richness",
    color = "Method",
    fill = "Method",
    shape = "Method"
  ) +
  theme(
    legend.title = element_blank()
  ) +
  scale_shape_manual(
    values = c(
      "CT_A" = 16,
      "CT_S" = 17,
      "PC"   = 15
    )
  ) +
  scale_color_manual(
    values = c(
      "CT_A" = "#CC79A7",
      "CT_S" = "#F0E442",
      "PC"   = "#666666"
    )
  ) +
  scale_fill_manual(
    values = c(
      "CT_A" = "#CC79A7",
      "CT_S" = "#F0E442",
      "PC"   = "#666666"
    )
  )

p_full

###################################################################
#extracting iNEXT estimates
#edit meta data to add location
# Rebuild PC only for method × site × location estimates

pc_detected_loc <- animals %>%
  filter(Method == "PC") %>%
  mutate(
    sample_id = paste("PC", Site, Location, Visit_No, PC_No, sep = "_")
  ) %>%
  distinct(Method, Site, Location, sample_id, Taxon)

pc_full_loc <- expand_grid(
  Site = unique(animals$Site[animals$Method == "PC"]),
  Location = unique(animals$Location[animals$Method == "PC"]),
  Visit_No = 1:3,
  PC_No = 1:10
) %>%
  mutate(
    Method = "PC",
    sample_id = paste("PC", Site, Location, Visit_No, PC_No, sep = "_")
  )

pc_all_loc <- pc_full_loc %>%
  left_join(
    pc_detected_loc,
    by = c("Method", "Site", "Location", "sample_id")
  )
#deal with CT
#rebuild CT detected to add location
ct_detected <- animals %>%
  filter(Method %in% c("CT_A", "CT_S")) %>%
  mutate(
    sample_id = paste(Method, Site, Camera_ID, Date, sep = "_")
  ) %>%
  dplyr::select(Method, Site, Camera_ID, Date, sample_id, Location, Taxon) %>%
  distinct()

#check
ct_all %>% count(Method, Location, is.na(Taxon))

#all
ct_all_loc <- ct_all %>%
  mutate(
    Group_ID = paste(Method, Site, Location, sep = "|")
  )

pc_all_loc <- pc_all_loc %>%
  mutate(
    Group_ID = paste(Method, Site, Location, sep = "|")
  )
# combine CT and PC
combined_loc <- bind_rows(
  ct_all_loc %>% dplyr::select(Group_ID, sample_id, Taxon),
  pc_all_loc %>% dplyr::select(Group_ID, sample_id, Taxon)
)

# make iNEXT input list
inext_group_list <- make_incidence_list(combined_loc)

# check groups
names(inext_group_list)
sapply(inext_group_list, function(x) x[1])
#iNEXT frequency-incidence input
make_incidence_list <- function(data) {
  
  effort <- data %>%
    distinct(Group_ID, sample_id) %>%
    count(Group_ID, name = "n_sampling_units")
  
  incidence <- data %>%
    filter(!is.na(Taxon), Taxon != "") %>%
    distinct(Group_ID, sample_id, Taxon) %>%
    group_by(Group_ID, Taxon) %>%
    summarise(incidence = n(), .groups = "drop")
  
  input <- incidence %>%
    left_join(effort, by = "Group_ID") %>%
    group_by(Group_ID) %>%
    summarise(
      data = list(c(first(n_sampling_units), incidence)),
      .groups = "drop"
    )
  
  setNames(input$data, input$Group_ID)
}


#run iNEXT
out_group <- iNEXT(
  inext_group_list,
  q = c(0, 1, 2),
  datatype = "incidence_freq"
)

estimated_method_site_location <- out_group$AsyEst %>%
  tidyr::separate(
    Assemblage,
    into = c("Method", "Site", "Location"),
    sep = "\\|",
    remove = FALSE
  )

writexl::write_xlsx(
  estimated_method_site_location,
  "estimated_iNEXT_method_site_location_nobirds.xlsx"
)

# iNEXT estimates by METHOD only

combined_method <- bind_rows(
  ct_all %>% dplyr::select(Method, sample_id, Taxon),
  pc_all %>% dplyr::select(Method, sample_id, Taxon)
) %>%
  mutate(Group_ID = Method)

method_list <- make_incidence_list(combined_method)

out_method <- iNEXT(
  method_list,
  q = c(0, 1, 2),
  datatype = "incidence_freq"
)

estimated_method_only <- out_method$AsyEst %>%
  rename(Method = Assemblage)

writexl::write_xlsx(
  estimated_method_only,
  "estimated_iNEXT_method_only_nobirds.xlsx"
)

estimated_method_only
############
###Rarefaction curves no birds

file_in  <- "overalluduplicatedtaxonomic_distinctness_input_template.xlsx"
file_out <- "overall_generated_units-birds.xlsx"

animals <- read_excel(file_in, sheet = "-birds") %>%
  mutate(
    Method = str_trim(as.character(Method)),
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT-S" ~ "CT_S",
      TRUE ~ Method
    ),
    Site = str_trim(as.character(Site)),
    Camera_ID = str_trim(as.character(Camera_ID)),
    Camera_ID = str_replace(Camera_ID, "cam$", ""),
    
    # FIX 2095 -> 2025
    Date = as.Date(`Day/Month`),
    Date = if_else(
      lubridate::year(Date) > 2050,
      Date - lubridate::years(70),
      Date
    ),
    
    Visit_No = as.integer(`Visit No.`),
    PC_No = as.integer(`PC No.`),
    Taxon = str_squish(paste(Genus, Species))
  ) %>%
  filter(
    Method %in% c("CT", "CT_A", "CT_S", "PC"),
    !is.na(Taxon),
    Taxon != "NA NA",
    Taxon != ""
  )

deploy <- read_excel(file_in, sheet = "camera_deployment_dates") %>%
  mutate(
    Method = str_trim(as.character(Method)),
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT-S" ~ "CT_S",
      TRUE ~ Method
    ),
    Site = str_trim(as.character(Site)),
    Camera_ID = str_trim(as.character(Camera_ID)),
    Camera_ID = str_replace(Camera_ID, "cam$", ""),
    start_date = as.Date(start_date),
    end_date = as.Date(end_date)
  )

# PC: generate all point-count units

pc_detected <- animals %>%
  filter(Method == "PC") %>%
  mutate(
    sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")
  ) %>%
  distinct(Method, Site, Visit_No, PC_No, sample_id, Taxon) %>%
  mutate(detected = 1)

pc_full <- expand_grid(
  Site = sort(unique(animals$Site[animals$Method == "PC"])),
  Visit_No = 1:3,
  PC_No = 1:10
) %>%
  mutate(
    Method = "PC",
    sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")
  )

pc_all <- pc_full %>%
  left_join(
    pc_detected,
    by = c("Method", "Site", "Visit_No", "PC_No", "sample_id")
  ) %>%
  mutate(
    detected = if_else(is.na(Taxon), 0L, 1L)
  ) %>%
  arrange(Site, Visit_No, PC_No, desc(detected), Taxon)


# CT: generate all active camera-day units

ct_full <- deploy %>%
  rowwise() %>%
  mutate(Date = list(seq(start_date, end_date, by = "day"))) %>%
  unnest(Date) %>%
  ungroup() %>%
  mutate(
    sample_id = paste(Method, Site, Camera_ID, Date, sep = "_")
  ) %>%
  select(Method, Site, Camera_ID, Date, sample_id)

ct_detected <- animals %>%
  filter(Method %in% c("CT", "CT_A", "CT_S")) %>%
  mutate(
    sample_id_detected = paste(Method, Site, Camera_ID, Date, sep = "_")
  ) %>%
  distinct(Method, Site, Camera_ID, Date, Taxon) %>%
  mutate(detected = 1)

ct_all <- ct_full %>%
  left_join(
    ct_detected %>%
      select(Site, Camera_ID, Date, Taxon, detected),
    by = c("Site", "Camera_ID", "Date")
  ) %>%
  mutate(
    detected = if_else(is.na(Taxon), 0L, 1L)
  ) %>%
  arrange(Method, Site, Camera_ID, Date, desc(detected), Taxon)

# Optional checks

pc_effort_check <- pc_all %>%
  distinct(Site, sample_id) %>%
  count(Site, name = "n_PC_units")

ct_effort_check <- ct_all %>%
  distinct(Method, Site, Camera_ID, sample_id) %>%
  count(Method, Site, Camera_ID, name = "n_CT_days")

zero_check <- bind_rows(
  pc_all %>% count(dataset = "PC", detected),
  ct_all %>% count(dataset = "CT", detected)
)


# Write generated sheets

write_xlsx(
  list(
    Generated_PC_units = pc_all,
    Generated_CT_units = ct_all,
    PC_effort_check = pc_effort_check,
    CT_effort_check = ct_effort_check,
    zero_detection_check = zero_check
  ),
  file_out
)

#rarefaction overall 

# Combine PC and CT generated sampling units
all_units <- bind_rows(
  pc_all %>% select(Method, Site, sample_id, Taxon, detected),
  ct_all %>% select(Method, Site, sample_id, Taxon, detected)
)

# Make one overall group
overall_units <- all_units %>%
  mutate(Group = "Overall")

# Build incidence-frequency input for iNEXT
overall_incidence <- overall_units %>%
  filter(!is.na(Taxon)) %>%
  distinct(Group, sample_id, Taxon) %>%
  group_by(Group, Taxon) %>%
  summarise(incidence = n(), .groups = "drop")

overall_effort <- overall_units %>%
  distinct(Group, sample_id) %>%
  count(Group, name = "n_sampling_units")

overall_input_df <- overall_incidence %>%
  left_join(overall_effort, by = "Group") %>%
  group_by(Group) %>%
  summarise(
    input = list(c(first(n_sampling_units), incidence)),
    .groups = "drop"
  )

overall_list <- setNames(
  overall_input_df$input,
  overall_input_df$Group
)

# Run iNEXT
overall_out <- iNEXT(
  overall_list,
  q = 0,
  datatype = "incidence_freq"
)

# Plot rarefaction curve
p_overall <- ggiNEXT(overall_out, type = 1) +
  theme_classic(base_size = 18) +
  labs(
    x = "Sampling units",
    y = "Species richness"
  ) +
  scale_color_manual(values = c("Overall" = "#2C7BB6")) +
  scale_fill_manual(values = c("Overall" = "#2C7BB6"))

p_overall


#RAREFACTION PER SITE

all_units <- bind_rows(
  pc_all %>% select(Method, Site, sample_id, Taxon, detected),
  ct_all %>% select(Method, Site, sample_id, Taxon, detected)
)

make_inext_input <- function(data, group_col) {
  
  data %>%
    filter(!is.na(Taxon)) %>%
    distinct({{ group_col }}, sample_id, Taxon) %>%
    count({{ group_col }}, Taxon, name = "incidence") %>%
    left_join(
      data %>%
        distinct({{ group_col }}, sample_id) %>%
        count({{ group_col }}, name = "n_units"),
      by = rlang::as_name(rlang::enquo(group_col))
    ) %>%
    group_by({{ group_col }}) %>%
    summarise(
      input = list(c(first(n_units), incidence)),
      .groups = "drop"
    ) %>%
    { setNames(.$input, .[[rlang::as_name(rlang::enquo(group_col))]]) }
}

site_list <- make_inext_input(all_units, Site)

site_out <- iNEXT(
  site_list,
  q = 0,
  datatype = "incidence_freq"
)

p_site <- ggiNEXT(site_out, type = 1, se = TRUE) +
  theme_classic(base_size = 18) +
  labs(
    x = "Number of sampling units",
    y = "Species richness",
    color = "Site",
    shape = "Site",
    linetype = "Curve type"
  ) +
  theme(
    legend.title = element_blank()
  ) +
  guides(
    fill= "none",
    color = guide_legend(
      order = 1,
      override.aes = list(
        linetype = 0,
        size = 3
      )
    ),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

p_site


#by location (shore vs inland)
pc_detected <- animals %>%
  filter(Method == "PC") %>%
  mutate(sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")) %>%
  distinct(Method, Site, Location, Visit_No, PC_No, sample_id, Taxon) %>%
  mutate(detected = 1)

pc_full <- animals %>%
  filter(Method == "PC") %>%
  distinct(Site, Location, Visit_No, PC_No) %>%
  mutate(
    Method = "PC",
    sample_id = paste("PC", Site, Visit_No, PC_No, sep = "_")
  )
pc_all <- pc_full %>%
  left_join(
    pc_detected,
    by = c("Method", "Site", "Location", "Visit_No", "PC_No", "sample_id")
  ) %>%
  mutate(
    detected = if_else(is.na(Taxon), 0L, 1L)
  )

ct_full <- deploy %>%
  rowwise() %>%
  mutate(Date = list(seq(start_date, end_date, by = "day"))) %>%
  unnest(Date) %>%
  ungroup() %>%
  mutate(sample_id = paste(Method, Site, Camera_ID, Date, sep = "_")) %>%
  select(Method, Site, Location, Camera_ID, Date, sample_id)
ct_detected <- animals %>%
  filter(Method %in% c("CT", "CT_A", "CT_S")) %>%
  distinct(Site, Camera_ID, Date, Taxon) %>%
  mutate(detected = 1)

ct_all <- ct_full %>%
  left_join(
    ct_detected,
    by = c("Site", "Camera_ID", "Date")
  ) %>%
  mutate(
    detected = if_else(is.na(Taxon), 0L, 1L)
  )
all_units_location_nobirds <- bind_rows(
  pc_all %>% select(sample_id, Taxon, Location),
  ct_all %>% select(sample_id, Taxon, Location)
) %>%
  filter(!is.na(Location))

location_list <- make_inext_input(all_units_location_nobirds, Location)

location_out <- iNEXT(
  location_list,
  q = 0,
  datatype = "incidence_freq"
)

p_location <- ggiNEXT(location_out, type = 1) +
  theme_classic(base_size = 18) +
  labs(
    x = "Number of sampling units",
    y = "Species richness",
    color = "Location",
    fill = "Location",
    shape = "Location",
    linetype = "Curve type"
  ) +
  theme(
    legend.title = element_blank()
  )

p_location

