library(readxl)
library(dplyr)
library(ggplot2)
library(DHARMa)
library(MASS)
library(ggeffects)
library(lme4)
library(glmmTMB)
library(tidyr)
library(lubridate)
library(readr)
library(car)
library(emmeans)
library(performance)


#GLM with offset effort
df <- read_excel("Exclusive records.xlsx")

# clean names 
names(df) <- c("Method", "Site", "Location", "Camera",
               "Records", "Effort", "Species", "Std_records")

df$Method <- factor(df$Method,
                    levels = c("CT_A", "CT_S"),
                    labels = c("CT_B", "CT_S"))
df$Site <- factor(df$Site)
df$Location <- factor(df$Location)

df$SiteLoc <- interaction(df$Site, df$Location)
#negative binomial

m_nb_glm <- glm.nb(
  Records ~ Method + offset(log(Effort)),
  data = df
)

summary(m_nb_glm)
#check assumptions
sim <- simulateResiduals(m_nb_glm)

plot(sim)

testDispersion(sim)

testZeroInflation(sim)

testUniformity(sim)
testOutliers(sim)

#omnibus
m_nb_null <- glm.nb(
  Records ~ offset(log(Effort)),
  data = df
)

anova(m_nb_null, m_nb_glm, test = "LRT")

#random effect model - tried, but fixed effects model was better

#lets plot

pred <- ggpredict(
  m_nb_glm,
  terms = "Method",
  condition = c(Effort = 10)
)
#combined
ggplot(df, aes(Method, Std_records, group = SiteLoc)) +
  geom_line(alpha = 0.3, color = "grey70") +
  geom_point(aes(color = Method), size = 2.5) +
  stat_summary(aes(group = Method),
               fun = mean,
               geom = "point",
               size = 5,
               color = "black") +
  scale_color_manual(
    values = c("CT_B" = "#CC79A7", "CT_S" = "#F0E442") 
    )+
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none"
  ) +
  labs(y = "Exclusive records per 10 trap days")

#Beta-binomial occupied vs empty

data <- read_excel("Cameracomparisonmodeling.xlsx",
                   sheet = "Sheet1"
)
data <- data %>%
  filter(!is.na(Site))
# Check column names
names(data)

# Convert grouping variables to factors
data$Camera <- factor(data$Camera)
data$Site <- factor(data$Site)

occ_data <- data %>%
  dplyr::transmute(
    Camera = factor(Camera),
    Site = factor(Site),
    
    occupied = as.numeric(`# of occupied videos`),
    empty = as.numeric(`Empty videos`)
  ) %>%
  dplyr::filter(
    !is.na(Camera),
    !is.na(Site),
    !is.na(occupied),
    !is.na(empty)
  )
nrow(occ_data)

str(occ_data)

sapply(occ_data, length)

summary(occ_data[, c("occupied", "empty")])

any(occ_data$occupied < 0)
any(occ_data$empty < 0)

occ_data %>%
  dplyr::mutate(total_videos = occupied + empty) %>%
  summary()

model_occ_bb <- glmmTMB(
  cbind(occupied, empty) ~ Camera + (1 | Site),
  family = betabinomial(link = "logit"),
  data = occ_data
)

summary(model_occ_bb)

#diagnositcs
set.seed(123)
sim_bb <- simulateResiduals(model_occ_bb)


plot(sim_bb)

testDispersion(sim_bb)

testUniformity(sim_bb)

testOutliers(sim_bb)

testZeroInflation(sim_bb)

#convert coefficents from log-odds to odds ratio
exp(fixef(model_occ_bb)$cond)


#get onmibus test
model_null <- glmmTMB(
  cbind(occupied, empty) ~ (1 | Site),
  family = betabinomial(link = "logit"),
  data = occ_data
)

anova(model_null, model_occ_bb)
#plot
occ_data <- data %>%
  transmute(
    Camera = factor(Camera),
    Site = factor(Site),
    Location = factor(Location),
    occupied = as.numeric(`# of occupied videos`),
    empty = as.numeric(`Empty videos`)
  ) %>%
  filter(
    Camera %in% c("CT_A", "CT_S"),
    !is.na(Site),
    !is.na(Location),
    !is.na(occupied),
    !is.na(empty),
    occupied + empty > 0
  ) %>%
  mutate(
    Camera = droplevels(Camera),
    SiteLoc = interaction(Site, Location),
    prop_occ = occupied / (occupied + empty)
  )

levels(occ_data$SiteLoc) <- c(
  "BWT Inland",
  "BWT Shore",
  "MVI Inland",
  "MVI Shore",
  "MVT Inland",
  "MVT Shore",
  "NFI Inland",
  "NFI Shore",
  "VLT Inland",
  "VLT Shore"
)

#cleveland by site and location
ggplot(
  occ_data,
  aes(
    x = prop_occ,
    y = SiteLoc,
    color = Camera
  )
) +
  geom_line(
    aes(group = SiteLoc),
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c(
      "CT_A" = "#CC79A7",
      "CT_S" = "#F0E442"
    ),
    labels = c(
      "CT_B",
      "CT_S"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic(base_size = 18) +
  labs(
    x = "Proportion of occupied videos",
    y = "Site and location",
    color = "Camera type"
  )


#cleveland by site
levels(occ_data$Site) <- c(
  "BWT",
  "MVI",
  "MVT",
  "NFI",
  "VLT"
)
ggplot(
  occ_data,
  aes(
    x = prop_occ,
    y = Site,
    color = Camera
  )
) +
  geom_line(
    aes(group = Site),
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c(
      "CT_A" = "#CC79A7",
      "CT_S" = "#F0E442"
    ),
    labels = c(
      "CT_B",
      "CT_S"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic(base_size = 18) +
  labs(
    x = "Proportion of occupied videos",
    y = "Site",
    color = "Camera type"
  )
#####
#bland-altman to compare between camera types

# Load data
detections <- read_excel("Combined_nohuman.xlsx", sheet = "Animals")
deploy <- read_excel("camera_deployment_dates.xlsx", sheet = "Sheet1")


# Clean detections

detections_clean <- detections %>%
  filter(Method %in% c("CT_N", "CT_S")) %>%
  mutate(
    date = as.Date(`Day/Month`),
    Camera_ID_clean = gsub("cam", "", Camera_ID),
    CameraType = case_when(
      Method == "CT_N" ~ "Battery",
      Method == "CT_S" ~ "Solar",
      TRUE ~ NA_character_
    ),
    Species_ID = paste(Genus, Species, sep = "_")
  )


# Clean deployment dates


deploy_clean <- deploy %>%
  mutate(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    Camera_ID_clean = Camera_ID,
    CameraType = case_when(
      Method %in% c("CT_A", "CT_N") ~ "Battery",
      Method == "CT_S" ~ "Solar",
      TRUE ~ NA_character_
    )
  )


# Get paired camera periods


paired_periods <- deploy_clean %>%
  select(Site, Location, CameraType, Camera_ID_clean, start_date, end_date) %>%
  pivot_wider(
    names_from = CameraType,
    values_from = c(Camera_ID_clean, start_date, end_date)
  ) %>%
  mutate(
    pair_start = pmax(start_date_Battery, start_date_Solar),
    pair_end = pmin(end_date_Battery, end_date_Solar)
  ) %>%
  filter(pair_start <= pair_end)


# Create one row per paired trap night


trap_nights <- paired_periods %>%
  rowwise() %>%
  mutate(date = list(seq(pair_start, pair_end, by = "day"))) %>%
  unnest(date) %>%
  ungroup() %>%
  select(
    Site, Location, date,
    Battery_ID = Camera_ID_clean_Battery,
    Solar_ID = Camera_ID_clean_Solar
  ) %>%
  pivot_longer(
    cols = c(Battery_ID, Solar_ID),
    names_to = "CameraType",
    values_to = "Camera_ID_clean"
  ) %>%
  mutate(
    CameraType = ifelse(CameraType == "Battery_ID", "Battery", "Solar")
  )


# Summarise detections per trap night


daily_detections <- detections_clean %>%
  group_by(Site, Location, date, Camera_ID_clean, CameraType) %>%
  summarise(
    detections = n(),
    richness = n_distinct(Species_ID),
    .groups = "drop"
  )

trapnight_summary <- trap_nights %>%
  left_join(
    daily_detections,
    by = c("Site", "Location", "date", "Camera_ID_clean", "CameraType")
  ) %>%
  mutate(
    detections = replace_na(detections, 0),
    richness = replace_na(richness, 0)
  )


# Bland-Altman data: detections

ba_detections <- trapnight_summary %>%
  select(Site, Location, date, CameraType, detections) %>%
  pivot_wider(
    names_from = CameraType,
    values_from = detections
  ) %>%
  mutate(
    mean_count = (Battery + Solar) / 2,
    difference = Battery - Solar
  )

bias_det <- mean(ba_detections$difference)
sd_det <- sd(ba_detections$difference)
loa_upper_det <- bias_det + 1.96 * sd_det
loa_lower_det <- bias_det - 1.96 * sd_det

bias_det
loa_upper_det
loa_lower_det

ggplot(ba_detections, aes(x = mean_count, y = difference)) +
  geom_count() +
  geom_hline(yintercept = bias_det, linewidth = 1) +
  geom_hline(yintercept = loa_upper_det, linetype = "dashed") +
  geom_hline(yintercept = loa_lower_det, linetype = "dashed") +
  labs(
    title = "Bland-Altman Plot: Detections per Trap Night",
    x = "Mean detections per trap night",
    y = "Difference: Battery - Solar"
  ) +
  theme_minimal(base_size = 18)


# Bland-Altman data: richness


ba_richness <- trapnight_summary %>%
  select(Site, Location, date, CameraType, richness) %>%
  pivot_wider(
    names_from = CameraType,
    values_from = richness
  ) %>%
  mutate(
    mean_richness = (Battery + Solar) / 2,
    difference = Battery - Solar
  )

bias_rich <- mean(ba_richness$difference)
sd_rich <- sd(ba_richness$difference)
loa_upper_rich <- bias_rich + 1.96 * sd_rich
loa_lower_rich <- bias_rich - 1.96 * sd_rich

bias_rich
loa_upper_rich
loa_lower_rich

ggplot(ba_richness, aes(x = mean_richness, y = difference)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = bias_rich, linewidth = 1) +
  geom_hline(yintercept = loa_upper_rich, linetype = "dashed") +
  geom_hline(yintercept = loa_lower_rich, linetype = "dashed") +
  labs(
    title = "Bland-Altman Plot: Species Richness per Trap Night",
    x = "Mean richness per trap night",
    y = "Difference: Battery - Solar"
  ) +
  theme_minimal(base_size = 18)

#final plot used to show how many were at same spot
ggplot(ba_richness, aes(x = mean_richness, y = difference)) +
  geom_count() +
  geom_hline(yintercept = bias_rich, linewidth = 1) +
  geom_hline(yintercept = loa_upper_rich, linetype = "dashed") +
  geom_hline(yintercept = loa_lower_rich, linetype = "dashed") +
  labs(
    title = "Bland–Altman Plot: Species Richness per Trap Night",
    x = "Mean richness per trap night",
    y = "Difference (Battery − Solar)"
  ) +
  theme_minimal(base_size = 18)

###
#quickly plot records per 10 SU from excel data

records <- data.frame(
  Method = factor(c("CT_B", "CT_S", "PC"),
                  levels = c("CT_B", "CT_S", "PC")),
  Records = c(7.1636, 10.6711, 12.8667)
)

ggplot(records, aes(x = Method, y = Records, fill = Method)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "CT_B" = "#CC79A7",
    "CT_S" = "#F0E442",
    "PC"   = "#666666"
  )) +
  labs(
    x = "Method",
    y = "Standardized records per 10 sampling units"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none"
  )

#################################
#cost method camera comparisons
#linear models Cost/record, TRIED Random effects in mixed models and log transforming, but linear model best

df <- read_excel(
  "Cameracomparisonmodeling.xlsx",
  sheet = "Sheet1"
) %>%
  rename(
    camera = Camera,
    site = Site,
    location = Location,
    cost_per_record = `Cost/record (EUR)`
  ) %>%
  mutate(
    camera = factor(camera),
    site = factor(site),
    location = factor(location),
    cost_per_record = parse_number(as.character(cost_per_record))
  )

df_cost <- df %>%
  filter(!is.na(cost_per_record))

#calculate descriptive stats
cost_summary <- df_cost %>%
  group_by(camera) %>%
  summarise(
    n = n(),
    mean = mean(cost_per_record),
    sd = sd(cost_per_record),
    se = sd/sqrt(n),
    median = median(cost_per_record),
    min = min(cost_per_record),
    max = max(cost_per_record)
  )

cost_summary

#box plot
ggplot(df_cost,
       aes(camera, cost_per_record, fill = camera)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  labs(
    x = "Camera type",
    y = "Cost per successful record (€)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

#linear model
m_cost <- lm(cost_per_record ~ camera,
             data = df_cost)

summary(m_cost)

#test assumptions
par(mfrow = c(2,2))
plot(m_cost)

#homogeneity of variance
car::leveneTest(cost_per_record ~ camera,
                data = df_cost)

emm_cost <- emmeans(m_cost, ~ camera)

pairs(
  emm_cost,
  adjust = "tukey"
)
############

