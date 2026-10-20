
# NMDS: Community composition by site, location, and method

library(readxl)
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(janitor)
library(dplyr)
library(ggplot2)
library(car)
library(DHARMa)
library(lmtest)
library(lme4)
library(permute)

# File path and analysis settings


file_path <- "methodduplicatedtaxonomic_distinctness_input_template.xlsx"


#use abundance
use_presence_absence <- FALSE

# Taxonomic resolution used in the NMDS.
taxonomic_level <- "Species"

# Random seed makes the NMDS reproducible.
set.seed(42)



# Import the two worksheets


animals_raw <- read_excel(
  path = file_path,
  sheet = "animals-domestic"
)

taxonomy_raw <- read_excel(
  path = file_path,
  sheet = "Taxonomy-domestic"
)


# Clean the animal detection records


animals <- animals_raw %>%
  clean_names() %>%
  transmute(
    Method      = str_trim(as.character(method)),
    Site        = str_trim(as.character(site)),
    Camera_ID   = str_trim(as.character(camera_id)),
    Date        = as.Date(day_month, origin = "1899-12-30"),
    Visit_No    = suppressWarnings(as.integer(visit_no)),
    PC_No       = suppressWarnings(as.integer(pc_no)),
    Location    = str_to_title(str_trim(as.character(location))),
    Common_name = str_squish(as.character(common_name)),
    Genus       = str_squish(as.character(genus)),
    Species     = str_squish(as.character(species))
  ) %>%
  mutate(
    # Construct the scientific name from genus and species.
    Taxon = str_squish(paste(Genus, Species)),
    
    # Standardise categorical values.
    Method = case_when(
      Method %in% c("CT_N", "CT-A", "CT_A") ~ "CT_A",
      Method %in% c("CT-S", "CT_S")         ~ "CT_S",
      Method == "PC"                        ~ "PC",
      TRUE                                  ~ Method
    ),
    
    Location = case_when(
      str_to_lower(Location) == "shore"  ~ "Shore",
      str_to_lower(Location) == "inland" ~ "Inland",
      TRUE                               ~ Location
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Method),
    !is.na(Location),
    !is.na(Genus),
    !is.na(Species),
    Genus != "",
    Species != ""
  )


# Check the categories found in the data.
animals %>%
  count(Site, Method, Location) %>%
  arrange(Site, Method, Location) %>%
  print(n = Inf)



# Clean and join the taxonomy worksheet


taxonomy <- taxonomy_raw %>%
  clean_names() %>%
  transmute(
    Taxon   = str_squish(as.character(taxon)),
    Species = str_squish(as.character(species)),
    Genus   = str_squish(as.character(genus)),
    Family  = str_squish(as.character(family)),
    Order   = str_squish(as.character(order)),
    Class   = str_squish(as.character(class)),
    Phylum  = str_squish(as.character(phylum)),
    Kingdom = str_squish(as.character(kingdom))
  ) %>%
  mutate(
    # Reconstructing Taxon helps eliminate extra spaces in the sheet.
    Taxon = str_squish(paste(Genus, Species))
  ) %>%
  distinct(Taxon, .keep_all = TRUE)

animals_tax <- animals %>%
  left_join(
    taxonomy %>%
      dplyr::select(
        Taxon,
        Taxonomy_Genus = Genus,
        Family,
        Order,
        Class,
        Phylum,
        Kingdom
      ),
    by = "Taxon"
  )


# Select taxonomic resolution


community_long <- animals_tax %>%
  mutate(
    Community_taxon = case_when(
      taxonomic_level == "Species" ~ Taxon,
      taxonomic_level == "Genus"   ~ Genus,
      taxonomic_level == "Family"  ~ Family,
      taxonomic_level == "Order"   ~ Order,
      taxonomic_level == "Class"   ~ Class,
      TRUE ~ Taxon
    )
  ) %>%
  filter(
    !is.na(Community_taxon),
    Community_taxon != ""
  )



# Create one sampling unit per Site × Method × Location


community_long <- community_long %>%
  mutate(
    Sample_ID = paste(Site, Method, Location, sep = "_")
  )


# Count the number of rows/detections for each taxon in each
# Site × Method × Location combination.
community_counts <- community_long %>%
  count(
    Sample_ID,
    Site,
    Method,
    Location,
    Community_taxon,
    name = "Abundance"
  )


# 7. Build the sample-by-taxon community matrix

community_wide <- community_counts %>%
  dplyr::select(
    Sample_ID,
    Community_taxon,
    Abundance
  ) %>%
  tidyr::pivot_wider(
    names_from = Community_taxon,
    values_from = Abundance,
    values_fill = 0
  )

# Convert the tibble into a numeric matrix.
community_matrix <- community_wide %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix()


# Store metadata separately in the same row order.
sample_metadata <- community_counts %>%
  distinct(Sample_ID, Site, Method, Location) %>%
  arrange(match(Sample_ID, rownames(community_matrix)))


stopifnot(
  identical(
    sample_metadata$Sample_ID,
    rownames(community_matrix)
  )
)


# Check the sampling units


cat("\nNumber of sampling units:", nrow(community_matrix), "\n")
cat("Number of taxa:", ncol(community_matrix), "\n")

sample_metadata %>%
  count(Site, Method, Location) %>%
  print(n = Inf)

# Sqrt transformation

community_matrix_analysis <- sqrt(community_matrix)

#environmental variables

# read data
model_data <- read_excel(
  "New model variables.xlsx",
  sheet = "-domestic"
) %>%
  clean_names()

names(model_data)

#vegetation PCA
# vegetation variables for PCA
veg_vars <- c(
  "floor_obstruction_percent",
  "understory_obstruction_percent",
  "nypa_density",      
  "presence_stilt_root"     
)

model_data_pca <- model_data %>%
  filter(if_all(all_of(veg_vars), ~ !is.na(.)))

pca_veg <- prcomp(
  model_data_pca[, veg_vars],
  center = TRUE,
  scale. = TRUE
)

summary(pca_veg)
pca_veg$rotation

model_data_pca$habitat_pc1 <- pca_veg$x[, 1]
model_data_pca$habitat_pc2 <- pca_veg$x[, 2]

#plot PCA
pca_scores <- as.data.frame(pca_veg$x)
pca_scores$location <- model_data_pca$location

pca_loadings <- as.data.frame(pca_veg$rotation)
pca_loadings$variable <- rownames(pca_loadings)

pca_loadings$label <- c(
  "Floor obstruction",
  "Understory obstruction",
  "Nypa abundance",
  "Prop root presence"
)

arrow_scale <- 2.5

ggplot(pca_scores, aes(PC1, PC2)) +
  
  geom_point(
    aes(color = location, shape = location),
    size = 4,
    alpha = 0.9
  ) +
  
  geom_segment(
    data = pca_loadings,
    aes(
      x = 0,
      y = 0,
      xend = PC1 * arrow_scale,
      yend = PC2 * arrow_scale
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.7
  ) +
  
  geom_text(
    data = pca_loadings,
    aes(
      x = PC1 * arrow_scale * 1.15,
      y = PC2 * arrow_scale * 1.15,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    nudge_y = c(0, 0, 0, -0.4)
  ) +
  
  scale_color_manual(values = c(
    "Inland" = "#D55E00",
    "Shore" = "#0072B2"
  )) +
  
  scale_shape_manual(values = c(
    "Inland" = 16,
    "Shore" = 17
  )) +
  
  theme_classic() +
  
  labs(
    x = paste0(
      "PC1 (",
      round(summary(pca_veg)$importance[2,1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_veg)$importance[2,2] * 100, 1),
      "%)"
    ),
    color = "Location",
    shape = "Location"
  )
#prepare model variables

model_data_pca <- model_data_pca %>%
  dplyr::mutate(
    site = factor(site),
    location = factor(location),
    
    habitat_pc1_z = as.numeric(scale(habitat_pc1)),
    habitat_pc2_z = as.numeric(scale(habitat_pc2)),
    
    tree_richness_z = as.numeric(scale(tree_species_richness)),
    ndvi_z = as.numeric(scale(ndvi)),
    
    human_activity_z = as.numeric(scale(human_activity)),
    shore_distance_z = as.numeric(scale(distance_from_shore_m)),
    dredging_distance_z = as.numeric(scale(distance_from_dredging_m))
  )

#standardize metadata names
community_metadata <- sample_metadata %>%
  dplyr::mutate(
    Method = dplyr::case_when(
      Method == "CT_A" ~ "CT_B",
      TRUE ~ as.character(Method)
    ),
    
    Site = dplyr::case_when(
      Site == "S1" ~ "DT",
      Site == "S2" ~ "DI",
      Site == "S3" ~ "BW",
      Site == "S4" ~ "VL",
      Site == "S5" ~ "NI",
      TRUE ~ as.character(Site)
    ),
    
    Location = stringr::str_to_title(
      as.character(Location)
    ),
    
    Sample_ID_clean = paste(
      Site,
      Method,
      Location,
      sep = "_"
    )
  )
rownames(community_matrix_analysis) <-
  community_metadata$Sample_ID_clean

stopifnot(
  identical(
    rownames(community_matrix_analysis),
    community_metadata$Sample_ID_clean
  )
)

#checking for matching identifiers in environmental data
unique(model_data_pca$site)
unique(model_data_pca$method)
unique(model_data_pca$location)

#standardize identifiers
env_data <- model_data_pca %>%
  dplyr::mutate(
    Site = dplyr::case_when(
      as.character(site) == "S1" ~ "DT",
      as.character(site) == "S2" ~ "DI",
      as.character(site) == "S3" ~ "BW",
      as.character(site) == "S4" ~ "VL",
      as.character(site) == "S5" ~ "NI",
      TRUE ~ as.character(site)
    ),
    
    Method = dplyr::case_when(
      method %in% c("CT_A", "CT-A", "CT_B") ~ "CT_B",
      method %in% c("CT_S", "CT-S") ~ "CT_S",
      method == "PC" ~ "PC",
      TRUE ~ as.character(method)
    ),
    
    Location = stringr::str_to_title(
      as.character(location)
    ),
    
    Sample_ID_clean = paste(
      Site,
      Method,
      Location,
      sep = "_"
    )
  )

#join environmental data to metadata
dbrda_metadata <- community_metadata %>%
  dplyr::select(
    Sample_ID_clean,
    Site,
    Method,
    Location
  ) %>%
  dplyr::left_join(
    env_data %>%
      dplyr::select(
        Sample_ID_clean,
        habitat_pc1_z,
        habitat_pc2_z,
        tree_richness_z,
        ndvi_z,
        human_activity_z,
        shore_distance_z,
        dredging_distance_z
      ),
    by = "Sample_ID_clean"
  )

#put metadata in same order
dbrda_metadata <- dbrda_metadata %>%
  dplyr::slice(
    match(
      rownames(community_matrix_analysis),
      Sample_ID_clean
    )
  )

#making one dataset
complete_rows <- stats::complete.cases(
  dbrda_metadata %>%
    dplyr::select(
      Method,
      Site,
      habitat_pc1_z,
      tree_richness_z,
      ndvi_z,
      human_activity_z,
      shore_distance_z,
      dredging_distance_z
    )
)

comm_dbrda <- community_matrix_analysis[
  complete_rows,
  ,
  drop = FALSE
]

env_dbrda <- dbrda_metadata[
  complete_rows,
  ,
  drop = FALSE
]

#remove units with no taxa
nonempty_rows <- rowSums(comm_dbrda) > 0

comm_dbrda <- comm_dbrda[
  nonempty_rows,
  ,
  drop = FALSE
]

env_dbrda <- env_dbrda[
  nonempty_rows,
  ,
  drop = FALSE
]

#fit canidate models
#habitat
dbrda_habitat <- vegan::dbrda(
  comm_dbrda ~
    habitat_pc1_z +
    tree_richness_z +
    ndvi_z +
    Condition(Method),
  data = env_dbrda,
  distance = "bray"
)

#disturbance
dbrda_disturbance <- vegan::dbrda(
  comm_dbrda ~
    human_activity_z +
    dredging_distance_z +
    Condition(Method),
  data = env_dbrda,
  distance = "bray"
)

#spatial
dbrda_spatial <- vegan::dbrda(
  comm_dbrda ~
    shore_distance_z +
    Condition(Method),
  data = env_dbrda,
  distance = "bray"
)

#site blocked overall permutation tests

perm_site <- how(
  nperm = 9999,
  blocks = env_dbrda$Site
)
#habiat
set.seed(42)
anova(
  dbrda_habitat,
  permutations = perm_site
)

#disturbance
set.seed(42)
anova(
  dbrda_disturbance,
  permutations = perm_site
)

#spatial
set.seed(42)
anova(
  dbrda_spatial,
  permutations = perm_site
)

#calculate explained variation
vegan::RsquareAdj(dbrda_habitat)
vegan::RsquareAdj(dbrda_disturbance)
vegan::RsquareAdj(dbrda_spatial)

#create summary table
dbrda_model_summary <- tibble::tibble(
  Model = c(
    "Habitat",
    "Disturbance",
    "Spatial"
  ),
  
  R2 = c(
    vegan::RsquareAdj(dbrda_habitat)$r.squared,
    vegan::RsquareAdj(dbrda_disturbance)$r.squared,
    vegan::RsquareAdj(dbrda_spatial)$r.squared
  ),
  
  Adjusted_R2 = c(
    vegan::RsquareAdj(dbrda_habitat)$adj.r.squared,
    vegan::RsquareAdj(dbrda_disturbance)$adj.r.squared,
    vegan::RsquareAdj(dbrda_spatial)$adj.r.squared
  )
)

dbrda_model_summary


