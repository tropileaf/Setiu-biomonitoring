library(readxl)
library(dplyr)
library(janitor)
library(ggplot2)
library(car)
library(DHARMa)
library(lmtest)
library(lme4)
library(emmeans)
library(glmmTMB)

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
  mutate(
    site = factor(site),
    
    habitat_pc1_z = as.numeric(scale(habitat_pc1)),
    habitat_pc2_z = as.numeric(scale(habitat_pc2)),
    
    tree_richness_z = as.numeric(scale(tree_species_richness)),
    ndvi_z = as.numeric(scale(ndvi)),
    
    human_activity_z = as.numeric(scale(human_activity)),
    shore_distance_z = as.numeric(scale(distance_from_shore_m)),
    dredging_distance_z = as.numeric(scale(distance_from_dredging_m))
  )
#look at correlation between variables
cor_data_z <- model_data_pca %>%
  dplyr::select(
    habitat_pc1_z,
    habitat_pc2_z,
    tree_richness_z,
    ndvi_z,
    human_activity_z,
    shore_distance_z,
    dredging_distance_z
  )

round(
  cor(cor_data_z, use = "complete.obs"),
  2
)

#assess multicolinearity
m_global_check <- glmmTMB::glmmTMB(
  observed_richness ~
    habitat_pc1_z +
    habitat_pc2_z +
    tree_richness_z +
    ndvi_z +
    human_activity_z +
    shore_distance_z +
    dredging_distance_z +
    offset(log(effort)),
  family = glmmTMB::nbinom2,
  data = model_data_pca
)

performance::check_collinearity(m_global_check)
#Negative binomial is best choice due to overdispersion in poisson 

#canidate models 
m_habitat <- glmmTMB(
  observed_richness ~
    habitat_pc1_z +
    tree_richness_z +
    ndvi_z +
    offset(log(effort)),
  family = nbinom2,
  data = model_data_pca
)

summary(m_habitat)

m_disturbance <- glmmTMB(
  observed_richness ~
    human_activity_z +
    dredging_distance_z +
    offset(log(effort)),
  family = nbinom2,
  data = model_data_pca
)

summary(m_disturbance)


m_spatial <- glmmTMB(
  observed_richness ~
    shore_distance_z +
    offset(log(effort)),
  family = nbinom2,
  data = model_data_pca
)

summary(m_spatial)

m_null <- glmmTMB(
  observed_richness ~
    offset(log(effort)),
  family = nbinom2,
  data = model_data_pca
)

summary(m_null)

AIC(m_habitat, m_disturbance, m_null, m_spatial)
#omnibus tests for canidate models
anova(m_null, m_habitat)
anova(m_null, m_disturbance)
anova(m_null, m_spatial)

#tried mixed model, but data doesn't support it
####check assumptions

#automate dharma plots
check_dharma <- function(model, model_name, n_sim = 1000) {
  
  set.seed(42)
  
  sim <- DHARMa::simulateResiduals(
    fittedModel = model,
    n = n_sim,
    plot = FALSE
  )
  
  cat("\n==============================\n")
  cat(model_name, "\n")
  cat("==============================\n")
  
  print(DHARMa::testUniformity(sim))
  print(DHARMa::testDispersion(sim))
  print(DHARMa::testZeroInflation(sim))
  print(DHARMa::testOutliers(sim))
  
  plot(sim)
  
  invisible(sim)
}

sim_habitat_nb <- check_dharma(
  m_habitat,
  "Including-birds habitat model"
)

sim_disturbance_nb <- check_dharma(
  m_disturbance,
  "Including-birds disturbance model"
)

sim_spatial_nb <- check_dharma(
  m_spatial,
  "Including-birds spatial model"
)

sim_null_nb <- check_dharma(
  m_null,
  "Including-birds null model"
)

#### Hill richness (q = 0): method comparisons ####

# Import and prepare data 
model_datamethod <- read_excel(
  "New model variables.xlsx",
  sheet = "-domestic"
) %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c("CT_A", "CT_S", "PC"),
      labels = c("CT_B", "CT_S", "PC")
    ),
    location = factor(location),
    site = factor(site),
    log_effort = log(effort)
  )


# Check variables used in the models
summary(model_datamethod[c(
  "observed_richness",
  "estimated_richness",
  "effort",
  "method"
)])

# Observed richness: Poisson GLM 

m_obs_pois <- glm(
  observed_richness ~ method + offset(log_effort),
  family = poisson(link = "log"),
  data = model_datamethod
)

summary(m_obs_pois)

# Overall effect of method
drop1(m_obs_pois, test = "Chisq")

# Estimated marginal means and Tukey comparisons
emm_obs <- emmeans(
  m_obs_pois,
  ~ method,
  type = "response"
)

summary(emm_obs)

pairs(
  emmeans(m_obs_pois, ~ method),
  adjust = "tukey"
)

# Diagnostics
set.seed(42)

sim_obs <- simulateResiduals(
  fittedModel = m_obs_pois,
  n = 1000,
  refit = FALSE,
  plot = FALSE
)

testResiduals(sim_obs)
testZeroInflation(sim_obs)
testQuantiles(sim_obs)

plot(sim_obs)


# Estimated richness: linear model 
m_est_method <- lm(
  log(estimated_richness + 1) ~ method,
  data = model_datamethod
)

summary(m_est_method)

# Estimated marginal means and Tukey comparisons
emm_est <- emmeans(
  m_est_method,
  ~ method
)

summary(emm_est)

pairs(
  emm_est,
  adjust = "tukey"
)

# Standard linear-model diagnostic plots
par(mfrow = c(2, 2))
plot(m_est_method)
par(mfrow = c(1, 1))


######

#prepare hill diversity data

library(tidyr)
hill_data <- read_excel(
  "estimated_iNEXT_method_site_location_corrected.xlsx",
  na = c("", "NA", "Inf", "-Inf")
) %>%
  clean_names() %>%
  mutate(
    observed = as.numeric(observed),
    estimator = as.numeric(estimator)
  )

str(hill_data[, c("observed", "estimator")])

hill_long <- hill_data %>%
  pivot_longer(
    cols = c(observed, estimator),
    names_to = "type",
    values_to = "diversity_value"
  ) %>%
  mutate(
    type = dplyr::recode(
      type,
      observed = "Observed diversity",
      estimator = "Estimated diversity"
    ),
    
    method = factor(
      method,
      levels = c("CT_A", "CT_S", "PC"),
      labels = c("CT_B", "CT_S", "PC")
    ),
    
    diversity_label = dplyr::recode(
      diversity,
      `Species richness` = "q = 0: Species richness",
      `Shannon diversity` = "q = 1: Shannon diversity",
      `Simpson diversity` = "q = 2: Simpson diversity"
    )
  )

hill_q0 <- hill_long %>%
  filter(diversity == "Species richness")

hill_q1 <- hill_long %>%
  filter(diversity == "Shannon diversity")

hill_q2 <- hill_long %>%
  filter(diversity == "Simpson diversity")

####plot
type_colors <- c(
  "Observed diversity"  = "#0072B2",
  "Estimated diversity" = "#D55E00"
)

ggplot(
  hill_long,
  aes(
    x = method,
    y = diversity_value,
    fill = type
  )
) +
  geom_boxplot(
    colour = "black",
    width = 0.78,
    alpha = 0.75,
    position = position_dodge(width = 0.82),
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_point(
    aes(colour = type),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.82
    ),
    size = 2.6,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ diversity_label,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = type_colors) +
  scale_colour_manual(values = type_colors) +
  scale_y_continuous(
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(
    ylim = c(0, 35),
    expand = FALSE
  ) +
  labs(
    x = "Method",
    y = "Effective number of species",
    fill = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(
      size = 18,
      margin = margin(t = 6)
    ),
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.7
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 18),
    legend.key.size = grid::unit(0.8, "cm"),
    legend.spacing.x = grid::unit(0.3, "cm"),
    panel.spacing.x = grid::unit(0.35, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha = 0.8,
        linewidth = 0.8,
        colour = "black"
      )
    ),
    colour = "none"
  )

# LM and pairwise for q1 
hill_q1_obs <- hill_data %>%
  filter(diversity == "Shannon diversity")

m_q1_obs <- lm(
  log(observed) ~ method,
  data = hill_q1_obs
)

#estimated
hill_q1_est <- hill_data %>%
  filter(diversity == "Shannon diversity")

m_q1_est <- lm(
  log(estimator) ~ method,
  data = hill_q1_est
)

# Pairwise comparisons and summary
summary(m_q1_obs)
emmeans(m_q1_obs, pairwise ~ method, adjust = "tukey")

summary(m_q1_est)
emmeans(m_q1_est, pairwise ~ method, adjust = "tukey")

#model and pairwise for simpson
hill_q2_obs <- hill_data %>%
  filter(diversity == "Simpson diversity")

m_q2_obs <- lm(
  log(observed) ~ method,
  data = hill_q2_obs
)

#estimated
hill_q2_est <- hill_data %>%
  filter(diversity == "Simpson diversity")

m_q2_est <- lm(
  log(estimator) ~ method,
  data = hill_q2_est
)

# Pairwise comparisons and summary
summary(m_q2_obs)
emmeans(m_q2_obs, pairwise ~ method, adjust = "tukey")

summary(m_q2_est)
emmeans(m_q2_est, pairwise ~ method, adjust = "tukey")

#assumptions
#shannon q=1
par(mfrow = c(2, 2))
plot(m_q1_obs)

par(mfrow = c(2, 2))
plot(m_q1_est)

#simpson q=2
par(mfrow = c(2, 2))
plot(m_q2_obs)

par(mfrow = c(2, 2))
plot(m_q2_est)

####random effect models - still singular fit
# Shannon observed
m_q1_obs_lmm <- lmer(
  log(observed) ~ method + (1 | site),
  data = hill_q1_obs,
  REML = TRUE
)

# Shannon estimated
m_q1_est_lmm <- lmer(
  log(estimator) ~ method + (1 | site),
  data = hill_q1_est,
  REML = TRUE
)

summary(m_q1_obs_lmm)
summary(m_q1_est_lmm)

####trying tweedie distribution - best was still log transformed LM
#shannon observed
hill_q1_obs <- hill_data %>%
  filter(diversity == "Shannon diversity")

m_q1_obs_tw <- glmmTMB(
  observed ~ method,
  family = tweedie(link = "log"),
  data = hill_q1_obs
)

summary(m_q1_obs_tw)

# Pairwise comparisons on original response scale
emmeans(
  m_q1_obs_tw,
  pairwise ~ method,
  type = "response",
  adjust = "tukey"
)

#shannon estimated
hill_q1_est <- hill_data %>%
  filter(diversity == "Shannon diversity")

m_q1_est_tw <- glmmTMB(
  estimator ~ method,
  family = tweedie(link = "log"),
  data = hill_q1_est
)

summary(m_q1_est_tw)

emmeans(
  m_q1_est_tw,
  pairwise ~ method,
  type = "response",
  adjust = "tukey"
)

#simpson observed
hill_q2_obs <- hill_data %>%
  filter(diversity == "Simpson diversity")

m_q2_obs_tw <- glmmTMB(
  observed ~ method,
  family = tweedie(link = "log"),
  data = hill_q2_obs
)

summary(m_q2_obs_tw)

emmeans(
  m_q2_obs_tw,
  pairwise ~ method,
  type = "response",
  adjust = "tukey"
)

#simpson estimated
hill_q2_est <- hill_data %>%
  filter(diversity == "Simpson diversity")

m_q2_est_tw <- glmmTMB(
  estimator ~ method,
  family = tweedie(link = "log"),
  data = hill_q2_est
)

summary(m_q2_est_tw)

emmeans(
  m_q2_est_tw,
  pairwise ~ method,
  type = "response",
  adjust = "tukey"
)

# Function to calculate adjusted AIC for log-transformed LM
adjusted_AIC <- function(model, original_values) {
  AIC(model) + 2 * sum(log(original_values))
}

# Shannon observed
adjusted_AIC(m_q1_obs, hill_q1_obs$observed)
AIC(m_q1_obs_tw)

# Shannon estimated
adjusted_AIC(m_q1_est, hill_q1_est$estimator)
AIC(m_q1_est_tw)

# Simpson observed
adjusted_AIC(m_q2_obs, hill_q2_obs$observed)
AIC(m_q2_obs_tw)

# Simpson estimated
adjusted_AIC(m_q2_est, hill_q2_est$estimator)
AIC(m_q2_est_tw)

comparison <- data.frame(
  Diversity = c(
    "Shannon observed",
    "Shannon estimated",
    "Simpson observed",
    "Simpson estimated"
  ),
  
  LM_adjusted_AIC = c(
    adjusted_AIC(m_q1_obs, hill_q1_obs$observed),
    adjusted_AIC(m_q1_est, hill_q1_est$estimator),
    adjusted_AIC(m_q2_obs, hill_q2_obs$observed),
    adjusted_AIC(m_q2_est, hill_q2_est$estimator)
  ),
  
  Tweedie_AIC = c(
    AIC(m_q1_obs_tw),
    AIC(m_q1_est_tw),
    AIC(m_q2_obs_tw),
    AIC(m_q2_est_tw)
  )
)

comparison$Delta_AIC <- (
  comparison$LM_adjusted_AIC -
    comparison$Tweedie_AIC
)

print(comparison)

#assumptions 
#shannon observed
res_q1_obs <- simulateResiduals(
  fittedModel = m_q1_obs_tw,
  n = 1000
)

plot(res_q1_obs)

testUniformity(res_q1_obs)
testDispersion(res_q1_obs)
testOutliers(res_q1_obs)

#shannon estimated
res_q1_est <- simulateResiduals(
  fittedModel = m_q1_est_tw,
  n = 1000
)

plot(res_q1_est)

testUniformity(res_q1_est)
testDispersion(res_q1_est)
testOutliers(res_q1_est)

# Simpson observed
res_q2_obs <- simulateResiduals(
  fittedModel = m_q2_obs_tw,
  n = 1000
)

plot(res_q2_obs)
testUniformity(res_q2_obs)
testDispersion(res_q2_obs)
testOutliers(res_q2_obs)

# Simpson estimated
res_q2_est <- simulateResiduals(
  fittedModel = m_q2_est_tw,
  n = 1000
)

plot(res_q2_est)
testUniformity(res_q2_est)
testDispersion(res_q2_est)
testOutliers(res_q2_est)

#fit gamma models
# Shannon observed
m_q1_obs_gamma <- glmmTMB(
  observed ~ method,
  family = Gamma(link = "log"),
  data = hill_q1_obs
)

# Shannon estimated
m_q1_est_gamma <- glmmTMB(
  estimator ~ method,
  family = Gamma(link = "log"),
  data = hill_q1_est
)

# Simpson observed
m_q2_obs_gamma <- glmmTMB(
  observed ~ method,
  family = Gamma(link = "log"),
  data = hill_q2_obs
)

# Simpson estimated
m_q2_est_gamma <- glmmTMB(
  estimator ~ method,
  family = Gamma(link = "log"),
  data = hill_q2_est
)

# Shannon observed
res_q1_obs_gamma <- simulateResiduals(
  m_q1_obs_gamma, n = 1000
)
plot(res_q1_obs_gamma)
testUniformity(res_q1_obs_gamma)
testDispersion(res_q1_obs_gamma)
testOutliers(res_q1_obs_gamma)

# Shannon estimated
res_q1_est_gamma <- simulateResiduals(
  m_q1_est_gamma, n = 1000
)
plot(res_q1_est_gamma)
testUniformity(res_q1_est_gamma)
testDispersion(res_q1_est_gamma)
testOutliers(res_q1_est_gamma)

# Simpson observed
res_q2_obs_gamma <- simulateResiduals(
  m_q2_obs_gamma, n = 1000
)
plot(res_q2_obs_gamma)
testUniformity(res_q2_obs_gamma)
testDispersion(res_q2_obs_gamma)
testOutliers(res_q2_obs_gamma)

# Simpson estimated
res_q2_est_gamma <- simulateResiduals(
  m_q2_est_gamma, n = 1000
)
plot(res_q2_est_gamma)
testUniformity(res_q2_est_gamma)
testDispersion(res_q2_est_gamma)
testOutliers(res_q2_est_gamma)

#fixing weird tweedie convergence
m_q2_obs_tw_bfgs <- update(
  m_q2_obs_tw,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS")
  )
)

m_q2_est_tw_bfgs <- update(
  m_q2_est_tw,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS")
  )
)

# Check convergence
m_q2_obs_tw_bfgs$fit$convergence
m_q2_obs_tw_bfgs$sdr$pdHess

m_q2_est_tw_bfgs$fit$convergence
m_q2_est_tw_bfgs$sdr$pdHess

#compare tweedie and gamma AIC
# Shannon observed
AIC(m_q1_obs_tw, m_q1_obs_gamma)

# Shannon estimated
AIC(m_q1_est_tw, m_q1_est_gamma)

# Simpson observed
AIC(m_q2_obs_tw, m_q2_obs_gamma)

# Simpson estimated
AIC(m_q2_est_tw, m_q2_est_gamma)

comparison <- data.frame(
  Response = c(
    "Shannon observed",
    "Shannon estimated",
    "Simpson observed",
    "Simpson estimated"
  ),
  
  Tweedie_AIC = c(
    AIC(m_q1_obs_tw),
    AIC(m_q1_est_tw),
    AIC(m_q2_obs_tw_bfgs),
    AIC(m_q2_est_tw_bfgs)
  ),
  
  Gamma_AIC = c(
    AIC(m_q1_obs_gamma),
    AIC(m_q1_est_gamma),
    AIC(m_q2_obs_gamma),
    AIC(m_q2_est_gamma)
  )
)

comparison$Delta_AIC <- (
  comparison$Tweedie_AIC -
    comparison$Gamma_AIC
)

print(comparison)


####
#### Hill richness (q = 0): method comparisons excluding birds 

#Import and prepare data 

model_datamethod_no_birds <- read_excel(
  "New model variables.xlsx",
  sheet = "-birds"
) %>%
  clean_names() %>%
  mutate(
    method = factor(
      method,
      levels = c("CT_A", "CT_S", "PC"),
      labels = c("CT_B", "CT_S", "PC")
    ),
    location = factor(location),
    site = factor(site),
    log_effort = log(effort)
  )

#Prepare richness figure 

combined_no_birds <- model_datamethod_no_birds %>%
  dplyr::select(
    method,
    observed_richness,
    estimated_richness
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      observed_richness,
      estimated_richness
    ),
    names_to = "type",
    values_to = "richness"
  ) %>%
  dplyr::mutate(
    type = dplyr::recode(
      type,
      observed_richness = "Observed richness",
      estimated_richness = "Estimated richness"
    ),
    type = factor(
      type,
      levels = c(
        "Observed richness",
        "Estimated richness"
      )
    )
  )


type_colors <- c(
  "Observed richness" = "#0072B2",
  "Estimated richness" = "#D55E00"
)


ggplot(
  combined_no_birds,
  aes(
    x = method,
    y = richness,
    fill = type
  )
) +
  geom_boxplot(
    alpha = 0.8,
    colour = "black",
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = type_colors) +
  theme_classic(base_size = 18) +
  labs(
    x = "Method",
    y = "Non-avian species richness",
    fill = "Richness type"
  )


#Observed richness: Poisson GLM 

m_obs_pois_no_birds <- glm(
  observed_richness ~ method + offset(log_effort),
  family = poisson(link = "log"),
  data = model_datamethod_no_birds
)

summary(m_obs_pois_no_birds)


# Overall effect of method
drop1(
  m_obs_pois_no_birds,
  test = "Chisq"
)

# Estimated marginal means on response scale
emm_obs_no_birds <- emmeans(
  m_obs_pois_no_birds,
  ~ method,
  type = "response"
)

summary(emm_obs_no_birds)


# Tukey-adjusted pairwise comparisons
pairs(
  emmeans(
    m_obs_pois_no_birds,
    ~ method
  ),
  adjust = "tukey"
)

# DHARMa diagnostics
set.seed(42)

sim_obs_no_birds <- simulateResiduals(
  fittedModel = m_obs_pois_no_birds,
  n = 1000,
  refit = FALSE,
  plot = FALSE
)

testResiduals(sim_obs_no_birds)
testZeroInflation(sim_obs_no_birds)
testQuantiles(sim_obs_no_birds)

plot(sim_obs_no_birds)

# Estimated richness: linear model

m_est_no_birds <- lm(
  log(estimated_richness + 1) ~ method,
  data = model_datamethod_no_birds
)

summary(m_est_no_birds)


# Standard linear-model diagnostic plots
par(mfrow = c(2, 2))
plot(m_est_no_birds)
par(mfrow = c(1, 1))


# Estimated marginal means
emm_est_no_birds <- emmeans(
  m_est_no_birds,
  ~ method
)

summary(emm_est_no_birds)


# Tukey-adjusted pairwise comparisons
pairs(
  emm_est_no_birds,
  adjust = "tukey"
)

#again tried random, but didn't need

####investigating ecological variables without birds

# read data
model_data <- read_excel(
  "New model variables.xlsx",
  sheet = "-birds"
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

#prepare model variables

model_data_pca <- model_data_pca %>%
  mutate(
    site = factor(site),
    
    habitat_pc1_z = as.numeric(scale(habitat_pc1)),
    habitat_pc2_z = as.numeric(scale(habitat_pc2)),
    
    tree_richness_z = as.numeric(scale(tree_species_richness)),
    ndvi_z = as.numeric(scale(ndvi)),
    
    human_activity_z = as.numeric(scale(human_activity)),
    shore_distance_z = as.numeric(scale(distance_from_shore_m)),
    dredging_distance_z = as.numeric(scale(distance_from_dredging_m))
  )

# Habitat model
m_habitat_pois <- glm(
  observed_richness ~
    habitat_pc1_z +
    tree_richness_z +
    ndvi_z +
    offset(log(effort)),
  family = poisson(link = "log"),
  data = model_data_pca
)

summary(m_habitat_pois)

# Disturbance model
m_disturbance_pois <- glm(
  observed_richness ~
    human_activity_z +
    dredging_distance_z +
    offset(log(effort)),
  family = poisson(link = "log"),
  data = model_data_pca
)

summary(m_disturbance_pois)

# Spatial model
m_spatial_pois <- glm(
  observed_richness ~
    shore_distance_z +
    offset(log(effort)),
  family = poisson(link = "log"),
  data = model_data_pca
)
summary(m_spatial_pois)

# Null model
m_null_pois <- glm(
  observed_richness ~
    offset(log(effort)),
  family = poisson(link = "log"),
  data = model_data_pca
)
summary(m_null_pois)
#omnibus tests
anova(m_null_pois, m_habitat_pois)
anova(m_null_pois, m_disturbance_pois)
anova(m_null_pois, m_spatial_pois)
#check assumptions with dharma
check_dharma <- function(model, model_name, n_sim = 1000) {
  
  set.seed(42)
  
  sim <- DHARMa::simulateResiduals(
    fittedModel = model,
    n = n_sim,
    plot = FALSE
  )
  
  cat("\n==============================\n")
  cat(model_name, "\n")
  cat("==============================\n")
  
  print(DHARMa::testUniformity(sim))
  print(DHARMa::testDispersion(sim))
  print(DHARMa::testZeroInflation(sim))
  print(DHARMa::testOutliers(sim))
  
  plot(sim)
  
  invisible(sim)
}

#run for each model
sim_habitat_pois <- check_dharma(
  m_habitat_pois,
  "Bird-excluded habitat model"
)

sim_disturbance_pois <- check_dharma(
  m_disturbance_pois,
  "Bird-excluded disturbance model"
)

sim_spatial_pois <- check_dharma(
  m_spatial_pois,
  "Bird-excluded spatial model"
)

sim_null_pois <- check_dharma(
  m_null_pois,
  "Bird-excluded null model"
)
