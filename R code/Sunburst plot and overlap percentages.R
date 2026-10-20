#all####
library(readxl)
library(dplyr)
library(tidyr)
library(metafor)
library(ggplot2)
library(viridis)
library(colorspace)
library(ggnewscale)
library(colorspace)
library(ggVennDiagram)

data_raw <- read_excel("Totalrecordsorder.xlsx",
                       sheet="all")

#clean
data <- data_raw %>%
  mutate(
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT_S" ~ "CT_S",
      Method == "PC" ~ "PC"
    )
  )
#Presence absence
data_pa <- data %>%
  group_by(Site, Order, Genus, Species, Method) %>%
  summarise(presence = 1, .groups = "drop") %>%
  pivot_wider(
    names_from = Method,
    values_from = presence,
    values_fill = 0
  )
#building dataset
make_pair <- function(df, method1, method2) {
  df %>%
    mutate(
      both = ifelse(.data[[method1]] == 1 & .data[[method2]] == 1, 1, 0),
      m1_only = ifelse(.data[[method1]] == 1 & .data[[method2]] == 0, 1, 0),
      m2_only = ifelse(.data[[method2]] == 1 & .data[[method1]] == 0, 1, 0),
      overall = 1,
      comparison = paste(method1, "vs", method2)
    )
}

pair_CA_CS <- make_pair(data_pa, "CT_A", "CT_S")
pair_CA_PC <- make_pair(data_pa, "CT_A", "PC")
pair_CS_PC <- make_pair(data_pa, "CT_S", "PC")

all_pairs <- bind_rows(pair_CA_CS, pair_CA_PC, pair_CS_PC)

#calculate LRR
datalrr <- all_pairs %>%
  mutate(
    a = overall - m1_only,
    b = overall - m2_only,
    c = m1_only + both,
    d = m2_only + both
  )

#effect size
effect <- escalc(
  measure = "RR",
  ai = a,
  bi = b,
  ci = c,
  di = d,
  data = datalrr
)

model <- rma(yi, vi, data = effect)
forest(model)

model_by_comp <- rma(yi, vi, mods = ~ comparison, data = effect)
forest(model_by_comp)

effect_order <- datalrr %>%
  group_by(Order, comparison) %>%
  summarise(
    a = sum(a),
    b = sum(b),
    c = sum(c),
    d = sum(d),
    .groups = "drop"
  )

effect_order_es <- escalc("RR", ai=a, bi=b, ci=c, di=d, data=effect_order)

model_order <- rma(yi, vi, data=effect_order_es)
forest(model_order)


summary_methods <- data_pa %>%
  summarise(
    CT_A = sum(CT_A),
    CT_S = sum(CT_S),
    PC = sum(PC),
    total = n()
  ) %>%
  pivot_longer(-total, names_to="method", values_to="count") %>%
  mutate(prop = count / total)


#trying sunburst

data_raw <- read_excel("Totalrecordsorder.xlsx",
                       sheet="all")

data_clean <- data_raw %>%
  mutate(
    Method = case_when(
      Method == "CT_N" ~ "CT-A",
      Method == "CT_S" ~ "CT-S",
      Method == "PC" ~ "PC",
      TRUE ~ Method
    ),
    taxon = paste(Genus, Species)
  )

# Presence/absence per taxon and order
data_pa <- data_clean %>%
  distinct(Order, taxon, Method) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Method,
    values_from = presence,
    values_fill = 0
  )

# Detection category
data_cat <- data_pa %>%
  mutate(
    n_methods = PC + `CT-A` + `CT-S`,
    category = case_when(
      PC == 1 & `CT-A` == 0 & `CT-S` == 0 ~ "PC only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 0 ~ "CT-A only",
      PC == 0 & `CT-A` == 0 & `CT-S` == 1 ~ "CT-S only",
      n_methods >= 2 ~ "Overlap"
    )
  )

# Counts by order and detection category
outer_counts <- data_cat %>%
  count(Order, category, name = "n")

# Order totals
order_totals <- data_cat %>%
  count(Order, name = "order_total")

# Overall totals
total_n <- nrow(data_cat)

# Center pie
overall <- data_cat %>%
  count(category, name = "n") %>%
  mutate(
    pct = n / sum(n) * 100,
    ymax = cumsum(n) / sum(n),
    ymin = lag(ymax, default = 0),
    label_y = (ymin + ymax) / 2
  )

# Inner order ring
inner <- order_totals %>%
  mutate(
    ymax = cumsum(order_total) / sum(order_total),
    ymin = lag(ymax, default = 0),
    label_y = (ymin + ymax) / 2
  )

orders <- unique(inner$Order)

order_colors <- setNames(
  viridis::viridis(length(orders)),
  orders
)

outer_colors <- setNames(
  colorspace::lighten(order_colors, amount = 0.45),
  orders
)
# Outer ring: segment widths still align with overall sunburst,
# but labels are now percentages within each order
outer <- outer_counts %>%
  left_join(order_totals, by = "Order") %>%
  left_join(inner %>% select(Order, order_ymin = ymin, order_ymax = ymax), by = "Order") %>%
  group_by(Order) %>%
  arrange(Order, category, .by_group = TRUE) %>%
  mutate(
    frac_within_order = n / order_total,
    ymin = order_ymin + cumsum(lag(frac_within_order, default = 0)) * (order_ymax - order_ymin),
    ymax = order_ymin + cumsum(frac_within_order) * (order_ymax - order_ymin),
    label_y = (ymin + ymax) / 2,
    pct = frac_within_order * 100,
    label = paste0(round(pct), "%"),
    method_symbol = case_when(
      category == "PC only" ~ "P",
      category == "CT-A only" ~ "A",
      category == "CT-S only" ~ "S",
      category == "Overlap" ~ "O"
    )
  ) %>%
  ungroup()

# Add method-specific shade variable
outer <- outer %>%
  mutate(
    outer_fill = paste(Order, category, sep = "_")
  )

# Create method-specific shades for each Order color
outer_colors_method <- unlist(lapply(names(order_colors), function(ord) {
  base_col <- order_colors[[ord]]
  
  cols <- c(
    "PC only"   = darken(base_col, 0.25),
    "CT-A only" = lighten(base_col, 0.10),
    "CT-S only" = lighten(base_col, 0.30),
    "Overlap"  = desaturate(lighten(base_col, 0.45), 0.3)
  )
  
  names(cols) <- paste(ord, names(cols), sep = "_")
  cols
}))

ggplot() +
  
  # center pie
  geom_rect(
    data = overall,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 1, fill = category),
    color = "white"
  ) +
  
  geom_text(
    data = overall,
    aes(x = 0.55, y = label_y, label = paste0(round(pct, 1), "%")),
    size = 3
  ) +
  
  scale_fill_manual(
    values = c(
      "PC only" = "grey20",
      "CT-A only" = "grey45",
      "CT-S only" = "grey75",
      "Overlap" = "grey92"
    ),
    name = "Method"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  
  # inner order ring
  geom_rect(
    data = inner,
    aes(ymin = ymin, ymax = ymax, xmin = 1.35, xmax = 2.55, fill = Order),
    color = "white"
  ) +
  
  geom_text(
    data = inner,
    aes(x = 1.95, y = label_y, label = Order),
    size = 3
  ) +
  
  scale_fill_manual(
    values = order_colors,
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # outer ring: lighter order colors
  
  geom_rect(
    data = outer,
    aes(ymin = ymin, ymax = ymax,
        xmin = 2.55, xmax = 3.25,
        fill = outer_fill),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = outer_colors_method,
    guide = "none"
  ) +
  # outer ring percentages
  geom_text(
    data = outer,
    aes(x = 3.43, y = label_y, label = label),
    size = 2.5
  ) +
  # method symbols on outside ring
  geom_point(
    data = outer,
    aes(x = 3.15, y = label_y, shape = category),
    size = 4,
    color = "black",
    fill = "white"
  ) +
  
  geom_text(
    data = outer,
    aes(x = 3.15, y = label_y, label = method_symbol),
    size = 2.3,
    fontface = "bold"
  ) +
  
  scale_shape_manual(
    values = c(
      "PC only" = 21,
      "CT-A only" = 21,
      "CT-S only" = 21,
      "Overlap" = 21
    ),
    labels = c(
      "PC only" = "Point count",
      "CT-A only" = "Affordable camera trap",
      "CT-S only" = "Solar camera trap",
      "Overlap" = "Overlap"
    ),
    name = "Detection category"
  ) +
  
  coord_polar(theta = "y") +
  xlim(0, 3.7) +
  theme_void()

####without any labels 

# Add method-specific shade variable
outer <- outer %>%
  mutate(
    outer_fill = paste(Order, category, sep = "_")
  )

# Create method-specific shades for each Order color
outer_colors_method <- unlist(lapply(names(order_colors), function(ord) {
  base_col <- order_colors[[ord]]
  
  cols <- c(
    "PC only"   = darken(base_col, 0.25),
    "CT-A only" = lighten(base_col, 0.10),
    "CT-S only" = lighten(base_col, 0.30),
    "Overlap"  = desaturate(lighten(base_col, 0.45), 0.3)
  )
  
  names(cols) <- paste(ord, names(cols), sep = "_")
  cols
}))
ggplot() +
  
  # center pie
  geom_rect(
    data = overall,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 1, fill = category),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = c(
      "PC only" = "grey20",
      "CT-A only" = "grey45",
      "CT-S only" = "grey75",
      "Overlap" = "grey92"
    ),
    name = "Method"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  
  # inner order ring
  geom_rect(
    data = inner,
    aes(ymin = ymin, ymax = ymax, xmin = 1.35, xmax = 2.55, fill = Order),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = order_colors,
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # outer ring: lighter order colors
  
  geom_rect(
    data = outer,
    aes(ymin = ymin, ymax = ymax,
        xmin = 2.55, xmax = 3.25,
        fill = outer_fill),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = outer_colors_method,
    guide = "none"
  ) +
  
  coord_polar(theta = "y") +
  xlim(0, 3.7) +
  theme_void()


#####trying venn and getting overlap percentages

# Create taxon lists for each method
venn_taxa <- list(
  PC = data_pa %>%
    filter(PC == 1) %>%
    pull(taxon),
  
  CT_B = data_pa %>%
    filter(`CT-A` == 1) %>%
    pull(taxon),
  
  CT_S = data_pa %>%
    filter(`CT-S` == 1) %>%
    pull(taxon)
)

# Plot Venn diagram
ggVennDiagram(
  venn_taxa,
  label = "count",
  label_alpha = 0,
  edge_size = 0.8,
  set_size = 5
) +
  scale_fill_gradient(low = "white", high = "grey50") +
  theme_void()


overlap_summary <- data_pa %>%
  mutate(
    overlap_group = case_when(
      PC == 1 & `CT-A` == 0 & `CT-S` == 0 ~ "PC only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 0 ~ "CT_B only",
      PC == 0 & `CT-A` == 0 & `CT-S` == 1 ~ "CT_S only",
      PC == 1 & `CT-A` == 1 & `CT-S` == 0 ~ "PC + CT_B only",
      PC == 1 & `CT-A` == 0 & `CT-S` == 1 ~ "PC + CT_S only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 1 ~ "CT_B + CT_S only",
      PC == 1 & `CT-A` == 1 & `CT-S` == 1 ~ "All three"
    )
  ) %>%
  count(overlap_group, name = "Taxa") %>%
  mutate(
    Percentage = 100 * Taxa / sum(Taxa)
  ) %>%
  arrange(desc(Taxa))

overlap_summary

#minus birds####
data_raw <- read_excel("Totalrecordsorder.xlsx",
                       sheet="-birds")

#clean
data <- data_raw %>%
  mutate(
    Method = case_when(
      Method == "CT_N" ~ "CT_A",
      Method == "CT_S" ~ "CT_S",
      Method == "PC" ~ "PC"
    )
  )
#Presence absence
data_pa <- data %>%
  group_by(Site, Order, Genus, Species, Method) %>%
  summarise(presence = 1, .groups = "drop") %>%
  pivot_wider(
    names_from = Method,
    values_from = presence,
    values_fill = 0
  )
#building dataset
make_pair <- function(df, method1, method2) {
  df %>%
    mutate(
      both = ifelse(.data[[method1]] == 1 & .data[[method2]] == 1, 1, 0),
      m1_only = ifelse(.data[[method1]] == 1 & .data[[method2]] == 0, 1, 0),
      m2_only = ifelse(.data[[method2]] == 1 & .data[[method1]] == 0, 1, 0),
      overall = 1,
      comparison = paste(method1, "vs", method2)
    )
}

pair_CA_CS <- make_pair(data_pa, "CT_A", "CT_S")
pair_CA_PC <- make_pair(data_pa, "CT_A", "PC")
pair_CS_PC <- make_pair(data_pa, "CT_S", "PC")

all_pairs <- bind_rows(pair_CA_CS, pair_CA_PC, pair_CS_PC)

#calculate LRR
datalrr <- all_pairs %>%
  mutate(
    a = overall - m1_only,
    b = overall - m2_only,
    c = m1_only + both,
    d = m2_only + both
  )

#effect size
effect <- escalc(
  measure = "RR",
  ai = a,
  bi = b,
  ci = c,
  di = d,
  data = datalrr
)

model <- rma(yi, vi, data = effect)
forest(model)

model_by_comp <- rma(yi, vi, mods = ~ comparison, data = effect)
forest(model_by_comp)

effect_order <- datalrr %>%
  group_by(Order, comparison) %>%
  summarise(
    a = sum(a),
    b = sum(b),
    c = sum(c),
    d = sum(d),
    .groups = "drop"
  )

effect_order_es <- escalc("RR", ai=a, bi=b, ci=c, di=d, data=effect_order)

model_order <- rma(yi, vi, data=effect_order_es)
forest(model_order)


summary_methods <- data_pa %>%
  summarise(
    CT_A = sum(CT_A),
    CT_S = sum(CT_S),
    PC = sum(PC),
    total = n()
  ) %>%
  pivot_longer(-total, names_to="method", values_to="count") %>%
  mutate(prop = count / total)


#trying sunburst

data_raw <- read_excel("Totalrecordsorder.xlsx",
                       sheet="-birds")

data_clean <- data_raw %>%
  mutate(
    Method = case_when(
      Method == "CT_N" ~ "CT-A",
      Method == "CT_S" ~ "CT-S",
      Method == "PC" ~ "PC",
      TRUE ~ Method
    ),
    taxon = paste(Genus, Species)
  )

# Presence/absence per taxon and order
data_pa <- data_clean %>%
  distinct(Order, taxon, Method) %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Method,
    values_from = presence,
    values_fill = 0
  )

# Detection category
data_cat <- data_pa %>%
  mutate(
    n_methods = PC + `CT-A` + `CT-S`,
    category = case_when(
      PC == 1 & `CT-A` == 0 & `CT-S` == 0 ~ "PC only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 0 ~ "CT-A only",
      PC == 0 & `CT-A` == 0 & `CT-S` == 1 ~ "CT-S only",
      n_methods >= 2 ~ "Overlap"
    )
  )

# Counts by order and detection category
outer_counts <- data_cat %>%
  count(Order, category, name = "n")

# Order totals
order_totals <- data_cat %>%
  count(Order, name = "order_total")

# Overall totals
total_n <- nrow(data_cat)

# Center pie
overall <- data_cat %>%
  count(category, name = "n") %>%
  mutate(
    pct = n / sum(n) * 100,
    ymax = cumsum(n) / sum(n),
    ymin = lag(ymax, default = 0),
    label_y = (ymin + ymax) / 2
  )

# Inner order ring
inner <- order_totals %>%
  mutate(
    ymax = cumsum(order_total) / sum(order_total),
    ymin = lag(ymax, default = 0),
    label_y = (ymin + ymax) / 2
  )

orders <- unique(inner$Order)

order_colors <- setNames(
  viridis::viridis(length(orders)),
  orders
)

outer_colors <- setNames(
  colorspace::lighten(order_colors, amount = 0.45),
  orders
)
# Outer ring: segment widths still align with overall sunburst,
# but labels are now percentages within each order
outer <- outer_counts %>%
  left_join(order_totals, by = "Order") %>%
  left_join(inner %>% select(Order, order_ymin = ymin, order_ymax = ymax), by = "Order") %>%
  group_by(Order) %>%
  arrange(Order, category, .by_group = TRUE) %>%
  mutate(
    frac_within_order = n / order_total,
    ymin = order_ymin + cumsum(lag(frac_within_order, default = 0)) * (order_ymax - order_ymin),
    ymax = order_ymin + cumsum(frac_within_order) * (order_ymax - order_ymin),
    label_y = (ymin + ymax) / 2,
    pct = frac_within_order * 100,
    label = paste0(round(pct), "%"),
    method_symbol = case_when(
      category == "PC only" ~ "P",
      category == "CT-A only" ~ "A",
      category == "CT-S only" ~ "S",
      category == "Overlap" ~ "O"
    )
  ) %>%
  ungroup()

# Add method-specific shade variable
outer <- outer %>%
  mutate(
    outer_fill = paste(Order, category, sep = "_")
  )

# Create method-specific shades for each Order color
outer_colors_method <- unlist(lapply(names(order_colors), function(ord) {
  base_col <- order_colors[[ord]]
  
  cols <- c(
    "PC only"   = darken(base_col, 0.25),
    "CT-A only" = lighten(base_col, 0.10),
    "CT-S only" = lighten(base_col, 0.30),
    "Overlap"  = desaturate(lighten(base_col, 0.45), 0.3)
  )
  
  names(cols) <- paste(ord, names(cols), sep = "_")
  cols
}))

ggplot() +
  
  # center pie
  geom_rect(
    data = overall,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 1, fill = category),
    color = "white"
  ) +
 
   geom_text(
    data = overall,
    aes(x = 0.55, y = label_y, label = paste0(round(pct, 1), "%")),
    size = 3
  ) +
 
   scale_fill_manual(
    values = c(
      "PC only" = "grey20",
      "CT-A only" = "grey45",
      "CT-S only" = "grey75",
      "Overlap" = "grey92"
    ),
    name = "Method"
  ) +
  
  ggnewscale::new_scale_fill() +
  

  # inner order ring
  geom_rect(
    data = inner,
    aes(ymin = ymin, ymax = ymax, xmin = 1.35, xmax = 2.55, fill = Order),
    color = "white"
  ) +
  
  geom_text(
    data = inner,
    aes(x = 1.95, y = label_y, label = Order),
    size = 3
  ) +
  
  scale_fill_manual(
    values = order_colors,
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # outer ring: lighter order colors
  
  geom_rect(
    data = outer,
    aes(ymin = ymin, ymax = ymax,
        xmin = 2.55, xmax = 3.25,
        fill = outer_fill),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = outer_colors_method,
    guide = "none"
  ) +
  # outer ring percentages
  geom_text(
    data = outer,
    aes(x = 3.43, y = label_y, label = label),
    size = 2.5
  ) +
  # method symbols on outside ring
  geom_point(
    data = outer,
    aes(x = 3.15, y = label_y, shape = category),
    size = 4,
    color = "black",
    fill = "white"
  ) +
  
  geom_text(
    data = outer,
    aes(x = 3.15, y = label_y, label = method_symbol),
    size = 2.3,
    fontface = "bold"
  ) +
  
  scale_shape_manual(
    values = c(
      "PC only" = 21,
      "CT-A only" = 21,
      "CT-S only" = 21,
      "Overlap" = 21
    ),
    labels = c(
      "PC only" = "Point count",
      "CT-A only" = "Affordable camera trap",
      "CT-S only" = "Solar camera trap",
      "Overlap" = "Overlap"
    ),
    name = "Detection category"
  ) +
  
  coord_polar(theta = "y") +
  xlim(0, 3.7) +
  theme_void()

####without any labels

# Add method-specific shade variable
outer <- outer %>%
  mutate(
    outer_fill = paste(Order, category, sep = "_")
  )

# Create method-specific shades for each Order color
outer_colors_method <- unlist(lapply(names(order_colors), function(ord) {
  base_col <- order_colors[[ord]]
  
  cols <- c(
    "PC only"   = darken(base_col, 0.25),
    "CT-A only" = lighten(base_col, 0.10),
    "CT-S only" = lighten(base_col, 0.30),
    "Overlap"  = desaturate(lighten(base_col, 0.45), 0.3)
  )
  
  names(cols) <- paste(ord, names(cols), sep = "_")
  cols
}))
ggplot() +
  
  # center pie
  geom_rect(
    data = overall,
    aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 1, fill = category),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = c(
      "PC only" = "grey20",
      "CT-A only" = "grey45",
      "CT-S only" = "grey75",
      "Overlap" = "grey92"
    ),
    name = "Method"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  
  # inner order ring
  geom_rect(
    data = inner,
    aes(ymin = ymin, ymax = ymax, xmin = 1.35, xmax = 2.55, fill = Order),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = order_colors,
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # outer ring: lighter order colors
  
  geom_rect(
    data = outer,
    aes(ymin = ymin, ymax = ymax,
        xmin = 2.55, xmax = 3.25,
        fill = outer_fill),
    color = "white"
  ) +
  
  scale_fill_manual(
    values = outer_colors_method,
    guide = "none"
  ) +

  coord_polar(theta = "y") +
  xlim(0, 3.7) +
  theme_void()


#####trying venn and getting overlap percentages

# Create taxon lists for each method
venn_taxa <- list(
  PC = data_pa %>%
    filter(PC == 1) %>%
    pull(taxon),
  
  CT_B = data_pa %>%
    filter(`CT-A` == 1) %>%
    pull(taxon),
  
  CT_S = data_pa %>%
    filter(`CT-S` == 1) %>%
    pull(taxon)
)

# Plot Venn diagram
ggVennDiagram(
  venn_taxa,
  label = "count",
  label_alpha = 0,
  edge_size = 0.8,
  set_size = 5
) +
  scale_fill_gradient(low = "white", high = "grey50") +
  theme_void()


overlap_summary <- data_pa %>%
  mutate(
    overlap_group = case_when(
      PC == 1 & `CT-A` == 0 & `CT-S` == 0 ~ "PC only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 0 ~ "CT_B only",
      PC == 0 & `CT-A` == 0 & `CT-S` == 1 ~ "CT_S only",
      PC == 1 & `CT-A` == 1 & `CT-S` == 0 ~ "PC + CT_B only",
      PC == 1 & `CT-A` == 0 & `CT-S` == 1 ~ "PC + CT_S only",
      PC == 0 & `CT-A` == 1 & `CT-S` == 1 ~ "CT_B + CT_S only",
      PC == 1 & `CT-A` == 1 & `CT-S` == 1 ~ "All three"
    )
  ) %>%
  count(overlap_group, name = "Taxa") %>%
  mutate(
    Percentage = 100 * Taxa / sum(Taxa)
  ) %>%
  arrange(desc(Taxa))

overlap_summary



