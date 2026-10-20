
# NMDS: community composition by method, site, and location

library(readxl)
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(janitor)
library(ggnewscale)
library(ggrepel)
library(indicspecies)

file_path <- "methodduplicatedtaxonomic_distinctness_input_template.xlsx"

#use abundance data
use_presence_absence <- FALSE

# Taxonomic resolution used in the NMDS.
taxonomic_level <- "Species"

# Random seed 
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

# Build the sample-by-taxon community matrix

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


# square-root transformation as Bray-Curtis can be strongly influenced by highly abundant taxa.

community_matrix_analysis <- sqrt(community_matrix)


# Run NMDS using Bray-Curtis dissimilarity

set.seed(42)

nmds_result <- metaMDS(
  community_matrix_analysis,
  distance      = "bray",
  k             = 2,
  trymax        = 999,
  autotransform = FALSE,
  wascores      = TRUE,
  trace         = TRUE
)

cat(
  "\nNMDS stress:",
  round(nmds_result$stress, 3),
  "\n"
)


# Extract NMDS site scores for plotting


nmds_scores <- scores(
  nmds_result,
  display = "sites",
  choices = c(1, 2)
) %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  as_tibble() %>%
  rename(
    NMDS1 = NMDS1,
    NMDS2 = NMDS2
  ) %>%
  left_join(sample_metadata, by = "Sample_ID")


# Format factors in the desired legend order.
nmds_scores <- nmds_scores %>%
  mutate(
    Site = factor(
      Site,
      levels = c("S1", "S2", "S3", "S4", "S5")
    ),
    Method = factor(
      Method,
      levels = c("PC", "CT_A", "CT_S")
    ),
    Location = factor(
      Location,
      levels = c("Shore", "Inland")
    ),
    Plot_label = paste(Site, Method, Location, sep = "_")
  )

#Faceted plot 
# Rename methods for plotting only
nmds_scores_plot <- nmds_scores %>%
  dplyr::mutate(
    
    Method = dplyr::case_when(
      Method == "CT_A" ~ "CT_B",
      TRUE ~ Method
    ),
    
    Site = factor(
      Site,
      levels = c("S1", "S2", "S3", "S4", "S5"),
      labels = c("MVT", "MVI", "BWT", "VLT", "NFI")
    )
  )

# Colours for methods
method_colours <- c(
  PC   = "#666666",
  CT_B = "#CC79A7",
  CT_S = "#F0E442"
)

# Shapes for sites
site_shapes <- c(
  MVT = 21,
  MVI = 22,
  BWT = 23,
  VLT = 24,
  NFI = 25
)

p_nmds_method <- ggplot(
  nmds_scores_plot,
  aes(x = NMDS1, y = NMDS2)
) +
  
  stat_ellipse(
    aes(
      group = Method,
      colour = Method
    ),
    level = 0.80,
    linetype = "dashed",
    linewidth = 0.8,
    type = "t",
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(
      fill = Method,
      shape = Site
    ),
    colour = "grey20",
    size = 5,
    stroke = 0.8,
    alpha = 0.9
  ) +
  
  facet_wrap(~Location) +
  
  scale_fill_manual(
    values = method_colours,
    name = "Method",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "grey20"
      )
    )
  ) +
  
  scale_colour_manual(
    values = method_colours,
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = site_shapes,
    name = "Site"
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.6,
    label = paste0(
      "Stress = ",
      round(nmds_result$stress, 3)
    ),
    size = 3.2,
    colour = "grey50",
    fontface = "italic"
  ) +
  
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS — assemblage composition coloured by method",
    subtitle = "Bray–Curtis | shape = site | panels = sampling location"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 18) +
  
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_nmds_method

####trying to get one simple figure
# Make sure the plotting variables are set up correctly
nmds_scores_plot <- nmds_scores %>%
  dplyr::mutate(
    Method = dplyr::case_when(
      Method == "CT_A" ~ "CT_B",
      TRUE ~ as.character(Method)
    ),
    
    Site = factor(
      Site,
      levels = c("S1", "S2", "S3", "S4", "S5"),
      labels = c("MVT", "MVI", "BWT", "VLT", "NFI")
    ),
    
    Method = factor(
      Method,
      levels = c("PC", "CT_B", "CT_S")
    ),
    
    Location = factor(
      Location,
      levels = c("Shore", "Inland")
    )
  )

site_colours <- c(
  MVT = "#0072B2",  # blue
  MVI = "#E69F00",  # orange
  BWT = "#CC79A7",  # reddish purple
  VLT = "#56B4E9",  # sky blue
  NFI = "#000000"   # black
)
method_shapes <- c(
  PC   = 21,
  CT_B = 22,
  CT_S = 24
)

#create shore and inland layers
p_all_clean <- ggplot2::ggplot(
  nmds_scores_plot,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2
  )
) +
  
  # Shore points: filled with the site colour
  ggplot2::geom_point(
    data = dplyr::filter(
      nmds_scores_plot,
      Location == "Shore"
    ),
    ggplot2::aes(
      shape = Method,
      colour = Site,
      fill = Site
    ),
    size = 5,
    stroke = 1
  ) +
  
  # Inland points: white fill with site-coloured borders
  ggplot2::geom_point(
    data = dplyr::filter(
      nmds_scores_plot,
      Location == "Inland"
    ),
    ggplot2::aes(
      shape = Method,
      colour = Site
    ),
    fill = "white",
    size = 5,
    stroke = 1
  ) +
  
  ggplot2::scale_colour_manual(
    values = site_colours,
    name = "Site"
  ) +
  
  ggplot2::scale_fill_manual(
    values = site_colours,
    guide = "none"
  ) +
  
  ggplot2::scale_shape_manual(
    values = method_shapes,
    name = "Method",
    labels = c(
      PC = "PC",
      CT_B = "CT_B",
      CT_S = "CT_S"
    )
  ) +
  
  ggplot2::annotate(
    geom = "text",
    x = Inf,
    y = -Inf,
    hjust = 1.08,
    vjust = -0.6,
    label = paste0(
      "Stress = ",
      round(nmds_result$stress, 3)
    ),
    size = 3.2,
    colour = "grey50",
    fontface = "italic"
  ) +
  
  ggplot2::labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS — assemblage composition",
    subtitle = paste0(
      "Bray–Curtis | colour = site | shape = method | ",
      "filled = shore, open = inland"
    )
  ) +
  
  ggplot2::coord_equal() +
  
  ggplot2::theme_classic(base_size = 13) +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

p_all_clean

location_legend_data <- tibble::tibble(
  NMDS1 = NA_real_,
  NMDS2 = NA_real_,
  Location = factor(
    c("Shore", "Inland"),
    levels = c("Shore", "Inland")
  )
)

p_all_clean <- p_all_clean +
  
  ggplot2::geom_point(
    data = location_legend_data,
    ggplot2::aes(
      x = NMDS1,
      y = NMDS2,
      alpha = Location
    ),
    shape = 21,
    size = 5,
    colour = "grey20",
    fill = c("grey20", "white"),
    inherit.aes = FALSE
  ) +
  
  ggplot2::scale_alpha_manual(
    values = c(
      Shore = 1,
      Inland = 1
    ),
    name = "Location",
    guide = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        fill = c("grey20", "white"),
        colour = "grey20",
        size = 5,
        alpha = 1
      )
    )
  )

p_all_clean

#trying with elipses for method
method_colours <- c(
  PC   = "#666666",
  CT_B = "#CC79A7",
  CT_S = "#F0E442"
)

p_all_method_ellipses <- ggplot2::ggplot(
  nmds_scores_plot,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2
  )
) +
  
  # Method ellipses
  ggplot2::stat_ellipse(
    ggplot2::aes(
      group = Method,
      colour = Method
    ),
    level = 0.80,
    type = "t",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.8,
    show.legend = TRUE
  ) +
  
  ggplot2::scale_colour_manual(
    values = method_colours,
    name = "Method ellipses"
  ) +
  
  # Start a second colour scale for the points
  ggnewscale::new_scale_colour() +
  
  # Shore: filled with site colour
  ggplot2::geom_point(
    data = dplyr::filter(
      nmds_scores_plot,
      Location == "Shore"
    ),
    ggplot2::aes(
      shape = Method,
      colour = Site,
      fill = Site
    ),
    size = 5,
    stroke = 1
  ) +
  
  # Inland: open with site-coloured outline
  ggplot2::geom_point(
    data = dplyr::filter(
      nmds_scores_plot,
      Location == "Inland"
    ),
    ggplot2::aes(
      shape = Method,
      colour = Site
    ),
    fill = "white",
    size = 5,
    stroke = 1
  ) +
  
  ggplot2::scale_colour_manual(
    values = site_colours,
    name = "Site"
  ) +
  
  ggplot2::scale_fill_manual(
    values = site_colours,
    guide = "none"
  ) +
  
  ggplot2::scale_shape_manual(
    values = method_shapes,
    name = "Method",
    labels = c(
      PC = "PC",
      CT_B = "CT-B",
      CT_S = "CT-S"
    )
  ) +
  # Dummy data for Location legend
  location_legend_data <- data.frame(
    NMDS1 = c(Inf, Inf),
    NMDS2 = c(Inf, Inf),
    Location = factor(c("Shore", "Inland"),
                      levels = c("Shore", "Inland"))
  )

p_all_method_ellipses <- p_all_method_ellipses +
  
  geom_point(
    data = location_legend_data,
    aes(
      x = NMDS1,
      y = NMDS2,
      alpha = Location
    ),
    shape = 21,
    size = 5,
    colour = "grey20",
    fill = c("grey20", "white"),
    inherit.aes = FALSE
  ) +
  
  scale_alpha_manual(
    values = c(
      Shore = 1,
      Inland = 1
    ),
    name = "Location",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill = c("grey20", "white"),
        colour = "grey20",
        size = 5,
        alpha = 1
      )
    )
  ) +
  ggplot2::annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.08,
    vjust = -0.6,
    label = paste0(
      "Stress = ",
      round(nmds_result$stress, 3)
    ),
    size = 3.2,
    colour = "grey50",
    fontface = "italic"
  ) +
  
  ggplot2::labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS — assemblage composition",
    subtitle = paste0(
      "Bray–Curtis | colour = site | shape = method | ",
      "filled = shore, open = inland | dashed ellipses = method"
    )
  ) +
  
  ggplot2::coord_equal() +
  
  ggplot2::theme_classic(base_size = 18) +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

p_all_method_ellipses


#EXTRACTING BRAY VALUES AND PERMANOVA
comm_for_tests <- sqrt(community_matrix)

sample_metadata_test <- sample_metadata %>%
  dplyr::mutate(
    Method = dplyr::case_when(
      as.character(Method) == "CT_A" ~ "CT_B",
      TRUE ~ as.character(Method)
    ),
    Method = factor(
      Method,
      levels = c("PC", "CT_B", "CT_S")
    ),
    Location = factor(
      Location,
      levels = c("Shore", "Inland")
    ),
    Site = factor(
      Site,
      levels = c("S1", "S2", "S3", "S4", "S5"),
      labels = c("DT", "DI", "BW", "VL", "NI")
    )
  )

stopifnot(
  identical(
    rownames(comm_for_tests),
    sample_metadata_test$Sample_ID
  )
)

#calculate disimilarity
bray_dist <- vegan::vegdist(
  comm_for_tests,
  method = "bray"
)

bray_dist

#blocked model to control for testing multiple times in sites
set.seed(42)

permanova_blocked <- vegan::adonis2(
  comm_for_tests ~ Method + Location,
  data = sample_metadata_test,
  method = "bray",
  permutations = 9999,
  by = "margin",
  strata = sample_metadata_test$Site
)

permanova_blocked

#test method site interaction
set.seed(42)

permanova_interaction <- vegan::adonis2(
  comm_for_tests ~ Method * Location,
  data = sample_metadata_test,
  method = "bray",
  permutations = 9999,
  by = "margin",
  strata = sample_metadata_test$Site
)

permanova_interaction

##pairwise for methods
pairwise_method_permanova <- function(
    community_matrix,
    metadata,
    permutations = 9999
) {
  
  method_pairs <- combn(
    levels(metadata$Method),
    2,
    simplify = FALSE
  )
  
  results <- lapply(
    method_pairs,
    function(method_pair) {
      
      keep <- metadata$Method %in% method_pair
      
      community_subset <- community_matrix[
        keep,
        ,
        drop = FALSE
      ]
      
      metadata_subset <- droplevels(
        metadata[keep, , drop = FALSE]
      )
      
      set.seed(42)
      
      test <- vegan::adonis2(
        community_subset ~ Method + Location,
        data = metadata_subset,
        method = "bray",
        permutations = permutations,
        by = "margin",
        strata = metadata_subset$Site
      )
      
      method_row <- test["Method", , drop = FALSE]
      
      tibble::tibble(
        Comparison = paste(
          method_pair,
          collapse = " vs "
        ),
        Df = method_row$Df,
        SumOfSqs = method_row$SumOfSqs,
        R2 = method_row$R2,
        F = method_row$F,
        p_value = method_row$`Pr(>F)`
      )
    }
  )
  
  dplyr::bind_rows(results) %>%
    dplyr::mutate(
      p_adjusted = stats::p.adjust(
        p_value,
        method = "BH"
      )
    )
}

pairwise_methods <- pairwise_method_permanova(
  community_matrix = comm_for_tests,
  metadata = sample_metadata_test,
  permutations = 9999
)

pairwise_methods

#shore inland (already done above)
shore_inland_result <- as.data.frame(
  permanova_blocked
) %>%
  tibble::rownames_to_column("Term") %>%
  dplyr::filter(Term == "Location")

shore_inland_result

#check dispersion
# Method dispersion
disp_method <- vegan::betadisper(
  bray_dist,
  group = sample_metadata_test$Method
)

set.seed(42)

method_dispersion_test <- vegan::permutest(
  disp_method,
  permutations = 9999
)

method_dispersion_test

# Location dispersion
disp_location <- vegan::betadisper(
  bray_dist,
  group = sample_metadata_test$Location
)

set.seed(42)

location_dispersion_test <- vegan::permutest(
  disp_location,
  permutations = 9999
)

location_dispersion_test

TukeyHSD(disp_method)

boxplot(
  disp_method,
  xlab = "Method",
  ylab = "Distance to centroid",
  main = "Multivariate dispersion among methods"
)
#export complete bray curtis matrix
bray_matrix <- as.matrix(bray_dist)

round(bray_matrix, 3)

write.csv(
  round(bray_matrix, 4),
  "Bray_Curtis_dissimilarity_matrix.csv",
  row.names = TRUE
)

#indicator species analysis

indval <- multipatt(
  comm_for_tests,
  sample_metadata_test$Method,
  control = how(nperm = 9999)
)

summary(indval)



