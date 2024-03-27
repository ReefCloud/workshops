## ---- packages
library(knitr)
library(tidyverse)
library(simstudy)
library(gstat)
library(sf)
library(RandomFields)
RFoptions(install = "no")
library(INLA)
library(inlabru)
library(mgcv)
library(stars)
library(patchwork)
library(Hmisc)
#library(downloadthis)
## ----end
 
## ---- set_seed
seed <- 123
## ----end

## Generate the spatial and temporal domains ==================================
## ---- temporal_domain
years <- 1:12
## ----end

## Create the full spatial domain
## ---- spatial_domain
spatial_domain <- st_geometry(
  st_multipoint(
    x = rbind(
      c(0, -10),
      c(3, -10),
      c(10, -20),
      c(1, -21),
      c(2, -16),
      c(0, -10)
    )
  )
) |>
  st_set_crs(4326) |>
  st_cast("POLYGON")

spatial_domain |>
  ggplot() +
  geom_sf() +
  theme_bw()
## ----end

## ---- spatial_grid
set.seed(seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(4236)

spatial_grid |> ggplot() +
  geom_sf(size = 0.1) +
  theme_bw()

## Compile the spatial data frame - note it is important for
## RFsimulate that the data be sorted by Longitude then Latitude
spatial_grid_pts_df <- spatial_grid |>
  st_coordinates() |>
  as.data.frame() |>
  dplyr::rename(Longitude = X, Latitude = Y) |>
  arrange(Longitude, Latitude)
## ----end

## Generate reefs ------------------------------------------------------
## Create a random field................................................
## ---- reefs_rf
RFoptions(seed = 1)
threshold <- 1.75
model <- RMexp(var = 1, scale = 0.1)
sim <- RFsimulate(model,
  x = as.vector(scale(spatial_grid_pts_df$Longitude,
    scale = FALSE
  )),
  y = as.vector(scale(spatial_grid_pts_df$Latitude,
    scale = FALSE
  ))
)
## combine with spatial data
reefs <- spatial_grid_pts_df |>
  mutate(Y = as.vector(sim))

reefs |>
  ggplot() +
  geom_tile(aes(y = Latitude, x = Longitude, fill = Y)) +
  coord_sf(crs = 4326) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )

reefs |>
  mutate(Y = Y > threshold) |>
  ggplot() +
  geom_tile(aes(y = Latitude, x = Longitude, fill = Y)) +
  coord_sf(crs = 4326) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
## ----end

## ---- reefs_rf_polygons
reefs_sf <- reefs |>
  st_as_sf(coords = c("Longitude", "Latitude")) |>
  filter(Y > threshold) |>
  st_buffer(0.05, endCapStyle = "SQUARE") |>
  st_cast("POLYGON") |>
  st_union() |>
  st_set_crs(4326)
reefs_sf |>
  ggplot() +
  geom_sf(data = spatial_domain, fill = NA) +
  geom_sf(fill = "red") + # nolint
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
reefs_full_sf <- reefs_sf
## ----end

## ---- reefs_rf_polygons_hollow
sf_use_s2(FALSE) # negative buffers dont work if this is true
reefs_sf <- reefs_sf |>
  st_buffer(0.01) |>
  st_difference(reefs_sf |> st_buffer(-0.01))
sf_use_s2(TRUE)
reefs_sf |> ggplot() +
  geom_sf() +
  coord_sf(xlim = c(2.4, 2.9), ylim = c(-16.75, -16.25)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
reefs_poly_sf <- reefs_sf |>
  st_cast("POLYGON") |>
  st_as_sf() |>
  mutate(Reef = paste0("Reef", 1:n()))
## ----end

## Generate synthetic (simulated) data ==================================
## The reponse data will be effected by the following:
## - base coral cover - the global average coral cover (pooled over space and time)
## - spatial pattern in this base cover which reflects the spatial pattern at T0
## - annual growth (e.g. 5-10% annual increase)
## - influence of covariates (spatio-temporal effects)
## - random noise

## ---- spatial_mesh
variance <- 1
kappa <- 1

alpha <- 2
mesh_pars <- c(1, 0.5, 0.1, 1, 0.5) *
  sqrt(alpha - ncol(spatial_grid_pts_df) / 2) / kappa
s <- inla.mesh.segment(
  spatial_grid_pts_df[chull(spatial_grid_pts_df), ]
)
mesh <- inla.mesh.2d(
  spatial_grid_pts_df[chull(spatial_grid_pts_df), ],
  max.edge = mesh_pars[1:2],
  cutoff = mesh_pars[3],
  offset = mesh_pars[4:5],
  boundary = s
)

ggplot() +
  gg(mesh) +
  geom_sf(data = spatial_domain, fill = NA, size = 2) +
  coord_sf(crs = 4326) +
  theme_bw()
## ----end

## ---- spatial_spde2
spde <- inla.spde2.matern(mesh, alpha = alpha)
## ----end

## ---- spatial_precision_matrix
theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
Q <- inla.spde2.precision(spde, theta = theta)
## ----end

## Calculate a lattice projection to and from the mesh
## ---- spatial_A
A <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(spatial_grid_pts_df)
)
# OR
## A <- inla.mesh.project(mesh = mesh,
##                        loc = as.matrix(spatial_grid_pts_df ))$A
## ----end

## Synthetic Baselines ----------------------------------------------------
## ---- baseline_spatial_HC
baseline_sample_hcc <- mesh$loc[, 1:2] |>
  as.data.frame() |>
  dplyr::select(Longitude = V1, Latitude = V2) |>
  mutate(
    clong = as.vector(scale(Longitude, scale = FALSE)),
    clat = as.vector(scale(Latitude, scale = FALSE)),
    Y = clong + sin(clat) + # rnorm(1,0,1) +
      1.5 * clong + clat
  ) |>
  mutate(Y = scales::rescale(Y, to = c(-2, 0.8)))

baseline_effects_hcc <- baseline_sample_hcc |>
  dplyr::select(Y) |>
  as.matrix()
baseline_pts_sample_hcc <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  baseline_effects_hcc
)
baseline_pts_effects_hcc <- baseline_pts_sample_hcc |>
  cbind() |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(Year = as.numeric(Year))

ggplot(baseline_pts_effects_hcc, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 7) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )

ggplot(baseline_pts_effects_hcc, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = 100 * plogis(Value))) +
  scale_fill_gradientn("Cover (%)", colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 7) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
## ----end

## ---- baseline_spatial_SC
baseline_sample_sc <- mesh$loc[, 1:2] |>
  as.data.frame() |>
  dplyr::select(Longitude = V1, Latitude = V2) |>
  mutate(
    clong = as.vector(scale(Longitude, scale = FALSE)),
    clat = as.vector(scale(Latitude, scale = FALSE)),
    Y = clong + sin(clat) + # rnorm(1,0,1) +
      1.5 * clong + -1.5 * clat
  ) |>
  mutate(Y = scales::rescale(Y, to = c(-4, -2)))

baseline_effects_sc <- baseline_sample_sc |>
  dplyr::select(Y) |>
  as.matrix()
baseline_pts_sample_sc <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  baseline_effects_sc
)
baseline_pts_effects_sc <- baseline_pts_sample_sc |>
  cbind() |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(Year = as.numeric(Year))

ggplot(baseline_pts_effects_sc, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 7) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
ggplot(baseline_pts_effects_sc, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = 100 * plogis(Value))) +
  scale_fill_gradientn("Cover (%)", colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 7) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
## ----end

## Synthetic DHW --------------------------------------------------------
## ---- covariate_DHW_temporal_trend
set.seed(seed)
dhw_temporal <- data.frame(Year = years) %>%
  mutate(
    cYear = Year - 1, # as.vector(scale(Year, scale=FALSE)),
    Y = 0.2 * cYear + sin(cYear),
    Y = Y * rbeta(length(years), Y, 1),
    Y = scales::rescale(Y - min(Y), to = c(0, 5))
  )
dhw_temporal %>%
  ggplot(aes(y = Y, x = Year)) +
  geom_line() +
  theme_bw(base_size = 7)
## ----end

## ---- covariate_DHW_effect
set.seed(seed)
dhw_sample <- inla.qsample(length(years),
  Q,
  seed = seed,
  constr = spde$f$extraconstr
)

rho <- rep(0.7, length(years))
rho <- rbeta(length(years), 0.2, 1)
x <- dhw_sample
for (j in 2:length(years)) {
  x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * dhw_sample[, j]
}
x <- sweep(x, 2, dhw_temporal$Y, FUN = "+")
dhw_effects <- scales::rescale(x, to = c(0, 1))
dhw_pts_sample <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  dhw_effects
)

dhw_pts_effects_df <- dhw_pts_sample %>%
  as.matrix() %>%
  as.data.frame() %>%
  cbind(spatial_grid_pts_df) %>%
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) %>%
  mutate(Year = as.numeric(Year))
## Value=scales::rescale(Value, to=c(0,1)))

## dhw.effect <- dhw.gmrf %>% mutate(Value=scales::rescale(Value, to=c(0,-0.3)))
ggplot(dhw_pts_effects_df, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, nrow = 2) +
  scale_fill_gradientn(colors = rev(heat.colors(10))) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())
## ----end

## Synthetic CYC --------------------------------------------------------
## ---- covariate_CYC_temporal_trend
set.seed(seed)
cyc <- vector("list", length(years))

## spatial_grid_pts_df |>
##     mutate(clong = as.vector(scale(Longitude, scale=FALSE)),
##            clat = as.vector(scale(Latitude, scale=FALSE)))

for (yr in years) {
  cat(paste("Year:", yr, "\n"))
  cyc_occur <- rbinom(1, 1, prob = min(0.05 * yr^2, 0.6))
  cat(paste("Cyclone Occurance:", cyc_occur, "\n"))
  cyc_intensity <- rbeta(1, 2, 1) |> round(2)
  cat(paste("Cyclone intensity:", cyc_intensity, "\n"))
  ## cyc_spatial <- spatial_grid_pts_df  |>
  lat_offset <- runif(1, 0, 5)
  cyc_spatial <- mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = lat_offset + runif(1, -1, 1) * clong + runif(1, -1, 1) *
        clat + sin(clat),
      # Y= Y - runif(1,-10,10),
      Y = abs(Y),
      Y = ifelse(Y > cyc_intensity, cyc_intensity, Y),
      Y = cyc_intensity - Y,
      Value = Y * cyc_occur
    )
  cyc[[yr]] <- cyc_spatial |> mutate(Year = yr)
}
cyc <- do.call("rbind", cyc)
cyc_effects_df <- cyc |>
  mutate(Value = scales::rescale(Value, to = c(0, 1)))

cyc_effects_df |>
  group_by(Year) |>
  summarise(
    Mean = mean(Value),
    Median = median(Value)
  ) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = "blue") +
  geom_line(aes(y = Median), color = "red") +
  theme_bw(base_size = 7)
## ----end

## ---- covariate_CYC_effect
cyc_effects <- cyc_effects_df |>
  dplyr::select(-clong, -clat, -Y) |>
  pivot_wider(
    id_cols = c(Longitude, Latitude),
    names_prefix = "sample:",
    names_from = Year,
    values_from = Value
  ) |>
  dplyr::select(-Longitude, -Latitude) |>
  as.matrix()

cyc_pts_sample <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  cyc_effects
)
cyc_pts_effects <- cyc_pts_sample |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(Year = as.numeric(Year))
## Value=scales::rescale(Value, to=c(0,-0.5)))

ggplot(cyc_pts_effects, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, nrow = 2) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
## ----end

## ---- CovariatesCYC_temporal_at_points
cyc_pts_effects |>
  group_by(Year) |>
  summarise(
    Mean = mean(Value),
    Median = median(Value)
  ) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = "blue") +
  geom_line(aes(y = Median), color = "red")
## ----end

## Synthetic Other --------------------------------------------------------
## ---- covariate_other_effect
set.seed(seed + 1)
other_sample <- inla.qsample(length(years),
  Q,
  seed = seed + 1,
  constr = spde$f$extraconstr
)

rho <- rep(0.7, length(years))
rho <- rbeta(length(years), 0.2, 1)
x <- other_sample
for (j in 2:length(years)) {
  x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * other_sample[, j]
}
other_effects <- scales::rescale(x, to = c(0, 1))
other_pts_sample <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  other_effects
)

other_pts_effects <- other_pts_sample |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(Year = as.numeric(Year)) # ,
## Value=scales::rescale(Value, to=c(0,1)))

ggplot(other_pts_effects, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, nrow = 2) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())
other_pts_effects |>
  group_by(Year) |>
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE)
  ) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean))
## ----end

## Compile all disturbances------------------------------------------------
## ---- compile_effects
disturb_effects <-
  (0.5 * dhw_effects) +
  (0.4 * cyc_effects) +
  (0.1 * other_effects) |>
  as.data.frame() # |>
all_effects_df <- mesh$loc[, 1:2] |>
  as.data.frame() |>
  dplyr::rename(Longitude = V1, Latitude = V2) |>
  cbind(disturb_effects) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = "Year",
    names_pattern = "sample:(.*)",
    values_to = "Y"
  ) |>
  mutate(Year = factor(Year, levels = sort(unique(as.numeric(
    as.character(Year)
  ))))) |>
  group_by(Longitude, Latitude) |>
  mutate(
    Growth_HCC = 0.3, ## Add growth onto this
    Growth_SC = 0.3,
    Y_HCC = cumsum(-Y + Growth_HCC), ## cumsum on link scale will accumulate effects
    Y_SC = cumsum(-Y + Growth_SC)
  )

## ----end

## ---- compile_effects_temporal
all_effects_df |>
  group_by(Year) |>
  summarise(
    Mean = mean(Y_HCC, na.rm = TRUE),
    Median = median(Y_HCC, na.rm = TRUE)
  ) |>
  ggplot(aes(x = as.numeric(as.character(Year)))) +
  geom_line(aes(y = Mean, color = "Mean")) +
  geom_line(aes(y = Median, color = "Median")) +
  scale_x_continuous("Year") +
  scale_y_continuous("Effect on HCC") +
  theme_bw()
## ----end

## ---- compile_effects_spatiotemporal
all_effects_df |>
  ggplot(aes(y = Latitude, x = Longitude)) +
  geom_point(aes(color = Y_HCC)) +
  ## geom_tile(aes(fill = Y_HCC)) +
  facet_wrap(~Year, nrow = 2) +
  scale_color_gradient2("HCC", low = "red", high = "green", mid = "white") +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())

all_effects <- all_effects_df |>
  pivot_wider(
    id_cols = c(Longitude, Latitude),
    names_prefix = "sample:",
    names_from = Year,
    values_from = Y_HCC
  )

## Project onto the spatial grid
disturb_pts_sample <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  all_effects |>
    ungroup() |> 
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()
)
disturb_pts_effects <- disturb_pts_sample |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(Year = as.numeric(Year))

ggplot(disturb_pts_effects, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, nrow = 2) +
  scale_fill_gradientn(colors = terrain.colors(10)) +
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1)
  )
## ----end

## Synthetic data ---------------------------------------------------------
## Hard coral
## ---- compile_synthetic_data_HCC
## Do all this on the link scale so that can use cumsum
all_effects_hcc <- all_effects_df |>
  full_join(baseline_sample_hcc |>
    dplyr::select(Longitude, Latitude, BASE_HCC = Y)) |>
  group_by(Longitude, Latitude) |>
  mutate(HCC = BASE_HCC + Y_HCC) |>
  ungroup() |>
  dplyr::select(-BASE_HCC, -Y_HCC) |>
  pivot_wider(
    id_cols = c(Longitude, Latitude),
    names_prefix = "sample:",
    names_from = Year,
    values_from = HCC
  ) |>
  dplyr::select(-Longitude, -Latitude) |>
  as.matrix()

all_pts_sample_hcc <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df[, 1:2]),
  all_effects_hcc
)
all_pts_effects_hcc <- all_pts_sample_hcc |>
  as.matrix() |>
  as.data.frame() |>
  cbind(spatial_grid_pts_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(
    Year = as.numeric(Year),
    ## Value = Value+qlogis(0.3))
    Value = Value
  )
## Value=scales::rescale(Value, to=c(-0.5,0)))
## ----end

## ---- compile_synthetic_data_HCC_temporal
all_pts_effects_hcc |>
  mutate(Value = plogis(Value)) |> 
  group_by(Year) |>
  summarise(Mean = mean(Value),
    Median = median(Value)) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = "blue") +
  geom_line(aes(y = Median), color = "red")
## ----end

## ---- compile_synthetic_data_HCC_spatiotemporal
all_pts_effects_hcc |>
  mutate(Value = plogis(Value)) |> 
  ggplot(aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  coord_sf(crs = 4236) +
  facet_wrap(~Year, nrow = 2) + 
  theme_bw(base_size = 7) +
  theme(axis.title = element_blank(),
    legend.position = c(0.95,0.95),
    legend.justification = c(1,1))
## ----end

## Soft coral
## ---- compile_synthetic_data_SC
## Do all this on the link scale so that can use cumsum
all_effects_sc <- all_effects_df |>
  full_join(baseline_sample_sc |> dplyr::select(Longitude, Latitude, BASE_SC=Y)) |>
  group_by(Longitude, Latitude) |>
  mutate(SC = BASE_SC + Y_SC) |>
  ungroup() |>
  dplyr::select(-BASE_SC, -Y_SC) |>
  pivot_wider(id_cols = c(Longitude, Latitude), 
    names_prefix = 'sample:',
    names_from = Year, 
    values_from = SC) |>
  dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

all_pts_sample_sc <- inla.mesh.project(mesh,
  loc = as.matrix(spatial_grid_pts_df [,1:2]),
  all_effects_sc)
all_pts_effects_sc = all_pts_sample_sc |> 
  as.matrix() |> 
  as.data.frame() |>
  cbind(spatial_grid_pts_df ) |> 
  pivot_longer(cols = c(-Longitude, -Latitude),
    names_to = c('Year'),
    names_pattern = 'sample:(.*)',
    values_to = 'Value') |>
  mutate(Year = as.numeric(Year),
    ## Value = Value+qlogis(0.3))
    Value = Value)
## Value=scales::rescale(Value, to=c(-0.5,0)))
## ----end

## ---- compile_synthetic_data_SC_temporal
all_pts_effects_sc |>
  mutate(Value = plogis(Value)) |> 
  group_by(Year) |>
  summarise(Mean = mean(Value),
    Median = median(Value)) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = "blue") +
  geom_line(aes(y = Median), color = "red")
## ----end

## ---- compile_synthetic_data_SC_spatiotemporal
all_pts_effects_sc |>
  mutate(Value = plogis(Value)) |> 
  ggplot(aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  coord_sf(crs = 4236) +
  facet_wrap(~Year, nrow = 2) + 
  theme_bw(base_size = 7) +
  theme(axis.title = element_blank(),
    legend.position = c(0.95,0.95),
    legend.justification = c(1,1))
## ----end

## Macroalgae
## ---- compile_synthetic_data_MA
## Do all this on the link scale so that can use cumsum
all_pts_effects_ma <- all_pts_effects_hcc |>
  dplyr::rename(HCC=Value) |> 
  bind_cols(all_pts_effects_sc |>
              dplyr::select(SC=Value)) |>
  mutate(Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
    ## MA = Total_Avail*rbeta(n(), 2, 1),
    MA = Total_Avail,
    Value = qlogis(MA)) |>
  dplyr::select(-HCC, -SC, -Total_Avail, -MA)
## ----end

## ---- compile_synthetic_data_MA_temporal
all_pts_effects_ma |>
  mutate(Value = plogis(Value)) |> 
  group_by(Year) |>
  summarise(Mean = mean(Value),
    Median = median(Value)) |>
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = "blue") +
  geom_line(aes(y = Median), color = "red")
## ----end

## ---- compile_synthetic_data_MA_spatiotemporal
all_pts_effects_ma |>
  mutate(Value = plogis(Value)) |> 
  ggplot(aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  coord_sf(crs = 4236) +
  facet_wrap(~Year, nrow = 2) + 
  theme_bw(base_size = 7) +
  theme(axis.title = element_blank(),
    legend.position = c(0.95,0.95),
    legend.justification = c(1,1))
## ----end

## Broad scale reef patterns ==================================================
## - rasterize the reefs frame
## - convert to points (centroids of raster cells)
## - filter to the values of 1
## - extract coordinates
## - convert to data frame

## ---- points_in_reefs
data_reefs_sf <- reefs_sf |>
  st_as_stars(dx = 0.01) |>  # rasterize
  st_as_sf(as_points = TRUE) |>
  filter(values == 1L)

data_reefs_df <- data_reefs_sf |>
  st_coordinates() |>
  as.data.frame() |>
  dplyr::rename(Longitude = X, Latitude = Y)
## ----end

## ---- project_onto_reefs_HCC
data_reefs_sample_hcc <- inla.mesh.project(mesh,
  loc = as.matrix(data_reefs_df[, 1:2]),
  all_effects_hcc
)
data_reefs_hcc <- data_reefs_sample_hcc |>
  as.matrix() |>
  as.data.frame() |>
  cbind(data_reefs_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(
    Year = as.numeric(Year),
    ## Value = Value + qlogis(0.3))
    Value = Value
  )

data_reefs_pts_hcc_sf <- data_reefs_hcc |>
  st_as_sf(coords = c("Longitude", "Latitude")) |>
  st_set_crs(st_crs(data_reefs_sf))
sf_use_s2(FALSE)
data_reefs_pts_hcc_sf <- data_reefs_pts_hcc_sf |>
  st_intersection(reefs_poly_sf)
sf_use_s2(TRUE)
## ----end

## ---- reefs_plot_HCC
sf_use_s2(FALSE)
data_reefs_pts_hcc_sf |>
  st_crop(xmin = 2.5, xmax = 3, ymin = -16.75, ymax = -16.25) |>
  mutate(Y = plogis(Value)) |>
  ggplot() +
  geom_sf(aes(color = Y)) +
  facet_wrap(~Year, nrow = 2) +
  scale_color_gradientn(colors = terrain.colors(10)) +
  coord_sf(
    crs = 4236,
    xlim = c(2.5, 3),
    ylim = c(-16.75, -16.25)
  ) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())
sf_use_s2(TRUE)
## ----end

## ---- project_onto_reefs_SC
data_reefs_sample_sc <- inla.mesh.project(mesh,
  loc = as.matrix(data_reefs_df[, 1:2]),
  all_effects_sc
)
data_reefs_sc <- data_reefs_sample_sc |>
  as.matrix() |>
  as.data.frame() |>
  cbind(data_reefs_df) |>
  pivot_longer(
    cols = c(-Longitude, -Latitude),
    names_to = c("Year"),
    names_pattern = "sample:(.*)",
    values_to = "Value"
  ) |>
  mutate(
    Year = as.numeric(Year),
    Value = Value
  )

data_reefs_pts_sc_sf <- data_reefs_sc |>
  st_as_sf(coords = c("Longitude", "Latitude")) |>
  st_set_crs(st_crs(data_reefs_sf))
sf_use_s2(FALSE)
data_reefs_pts_sc_sf <- data_reefs_pts_sc_sf |>
  st_intersection(reefs_poly_sf)
sf_use_s2(TRUE)
## ----end

## ---- reefs_plot_SC
sf_use_s2(FALSE)
data_reefs_pts_sc_sf |>
    st_crop(xmin = 2.5, xmax = 3, ymin = -16.75, ymax = -16.25) |>
    mutate(Y = plogis(Value)) |>
    ggplot() +
    geom_sf(aes(color = Y)) +
    facet_wrap(~Year, nrow = 2) +
    scale_color_gradientn(colors = terrain.colors(10)) + 
    coord_sf(crs = 4236, 
             xlim = c(2.5,3), 
             ylim = c(-16.75,-16.25)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_blank())
sf_use_s2(TRUE)
## ----end

## ---- project_onto_reefs_MA
data_reefs_ma <- data_reefs_hcc |>
  rename(HCC = Value) |>
  full_join(data_reefs_sc |> rename(SC = Value)) |>
  mutate(
    Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
    MA = Total_Avail,
    Value = qlogis(MA)
  ) |>
  dplyr::select(-HCC, -SC, -Total_Avail, -MA)

data_reefs_pts_ma_sf <- data_reefs_ma |>
  st_as_sf(coords = c("Longitude", "Latitude")) |>
  st_set_crs(st_crs(data_reefs_sf))
sf_use_s2(FALSE)
data_reefs_pts_ma_sf <- data_reefs_pts_ma_sf |>
  st_intersection(reefs_poly_sf)
sf_use_s2(TRUE)
## ----end

## ---- reefs_plot_MA
sf_use_s2(FALSE)
data_reefs_pts_ma_sf |>
  st_crop(xmin = 2.5, xmax = 3, ymin = -16.75, ymax = -16.25) |>
  mutate(Y = plogis(Value)) |>
  ggplot() +
  geom_sf(aes(color = Y)) +
  facet_wrap(~Year, nrow = 2) +
  scale_color_gradientn(colors = terrain.colors(10)) +
  coord_sf(
    crs = 4236,
    xlim = c(2.5, 3),
    ylim = c(-16.75, -16.25)
  ) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())
sf_use_s2(TRUE)
## ----end

## Sampling designs (large scale components) ==================================

## ---- project_onto_reefs
data_reefs_pts_sf <-
  data_reefs_pts_hcc_sf |>
  rename(HCC = Value) |>
  bind_cols(data_reefs_pts_sc_sf |>
              dplyr::select(SC = Value) |> st_drop_geometry()) |>
  bind_cols(data_reefs_pts_ma_sf |>
              dplyr::select(MA = Value) |> st_drop_geometry())
## ----end

## ---- number_of_locs
nLocs <- 25 # Number of 'reefs'
nSites <- 2 # Number of 'sites' within 'reefs'
## ----end

## INLA ------------------------------------------------------------

## Fixed locs ....................................................
## ---- fixed_locs
set.seed(seed)
## Start by randomly selecting nLocs Reefs
Reefs_fixed <- data_reefs_pts_sf |>
  st_drop_geometry() |>
  dplyr::select(Reef) |>
  distinct() |>
  sample_n(size = nLocs) |>
  pull(Reef)
## Then filter to these Reefs before selecting a single location within
## each of the Reefs
data_fixed_locs_sf <- data_reefs_pts_sf |>
  filter(Reef %in% Reefs_fixed) |>
  dplyr::select(Reef, geometry) |>
  distinct(.keep_all = TRUE) |> 
  group_by(Reef) |>
  sample_n(nSites) |>
  mutate(Site = paste0("S", 1:n())) |>
  ungroup() |> 
  st_join(data_reefs_pts_sf |> 
            dplyr::select(-Reef))
## ----end

## ---- fixed_locs_bubble
g <-
  data_fixed_locs_sf |>
  pivot_longer(cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value") |>
  group_by(Group) |>
  nest() |>  
  mutate(G = purrr::map2(.x=data, .y=Group,
    .f=function(.x, .y) {
      .x |> st_as_sf() |>
        mutate(Value = plogis(Value)) |>
        ggplot() +
        geom_sf(data = reefs_poly_sf, color = "lightgray") +
        geom_sf(data = spatial_domain, fill = NA) +
        geom_sf(aes(color = Value, size = Value)) +
        scale_color_gradientn(colors = terrain.colors(10)) + 
        facet_wrap(~Year, nrow = 2) +
        coord_sf(crs = 4236) +
        theme_bw(base_size = 12) +
        theme(axis.title = element_blank()) +
        ggtitle(.y)
    }))

g$G |>
  patchwork::wrap_plots(nrow = 3) 
## ----end

## ---- fixed_locs_means
data_reefs_sum <- data_reefs_pts_sf |>
  st_drop_geometry() |>
  group_by(Year) |>
  dplyr::summarize(across(c(HCC, SC, MA), mean)) |>
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Reef_mean"
  )
data_fixed_locs_sum <- data_fixed_locs_sf |>
  st_drop_geometry() |>
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value"
  ) |>
  group_by(Year, Group) |>
  dplyr::summarize(mean_cl_boot(Value))

g <-
  data_reefs_sum |>
  full_join(data_fixed_locs_sum) |>
  group_by(Group) |>
  nest() |>
  mutate(G = purrr::map2(
    .x = data, .y = Group,
    .f = function(.x, .y) {
      .x |>
        mutate(across(c(y, ymin, ymax, Reef_mean), plogis)) |>
        ggplot() +
        geom_line(aes(y = Reef_mean, x = Year, color = "Reef level")) +
        geom_point(aes(y = Reef_mean, x = Year, color = "Reef level")) +
        geom_line(aes(y = y, x = Year, color = "Site level")) +
        geom_ribbon(aes(
          y = y, x = Year, ymin = ymin, ymax = ymax,
          fill = "Site level"
        ), alpha = 0.2) +
        scale_y_continuous("Cover", labels = function(x) x * 100) +
        scale_color_manual("",
          breaks = c("Reef level", "Site level"),
          values = c("red", "blue"), limits = c("Reef level", "Site level")
        ) +
        scale_fill_manual("", breaks = c("Reef level", "Site level"), values = c("red", "blue"), limits = c("Reef level", "Site level")) +
        theme_bw() +
        ggtitle(.y)
    }
  ))

g$G |> patchwork::wrap_plots(nrow = 3)
## ----end

## Random locs ....................................................
## ---- random_locs
set.seed(seed)
## Start by randomly selecting nLocs different Reefs per year
Reefs_random <- data_reefs_pts_sf |>
  st_drop_geometry() |>
  dplyr::select(Reef, Year) |>
  distinct() |>
  group_by(Year) |>
  sample_n(nLocs) |>
  ungroup()
## Then select a single location within each of the reefs per year
data_random_locs_sf <- data_reefs_pts_sf |>
  right_join(Reefs_random) |>
  group_by(Reef, Year) |>
  sample_n(nSites) |>
  mutate(Site = paste0("S", 1:n())) |>
  arrange(Year, Reef) |>
  ungroup()
## ----end

## ---- random_locs_bubble
g <-
  data_random_locs_sf |>
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value"
  ) |>
  group_by(Group) |>
  nest() |>
  mutate(G = purrr::map2(.x = data, .y = Group,
    .f = function(.x, .y) {
      .x |>
        st_as_sf() |>
        mutate(Value = plogis(Value)) |>
        ggplot() +
        geom_sf(data = reefs_poly_sf, color = "lightgray") +
        geom_sf(data = spatial_domain, fill = NA) +
        geom_sf(aes(color = Value, size = Value)) +
        scale_color_gradientn(colors = terrain.colors(10)) +
        facet_wrap(~Year, nrow = 2) +
        coord_sf(crs = 4236) +
        theme_bw(base_size = 12) +
        theme(axis.title = element_blank()) +
        ggtitle(.y)
    }))

g$G |>
  patchwork::wrap_plots(nrow = 3)
## ----end

## ---- random_locs_means
data_reefs_sum <- data_reefs_pts_sf |>
  st_drop_geometry() |>
  group_by(Year) |>
  summarize(across(c(HCC, SC, MA), mean)) |>
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Reef_mean"
  )

data_random_locs_sum <- data_random_locs_sf |>
  st_drop_geometry() |>
  pivot_longer(cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value") |>
  group_by(Year, Group) |>
  summarize(mean_cl_boot(Value))

g <-
  data_reefs_sum |>
  full_join(data_random_locs_sum) |>
  group_by(Group) |>
  nest() |>
  mutate(G = purrr::map2(
    .x = data, .y = Group,
    .f = function(.x, .y) {
      .x |>
        mutate(across(c(y, ymin, ymax, Reef_mean), plogis)) |>
        ggplot() +
        geom_line(aes(y = Reef_mean, x = Year, color = "Reef level")) +
        geom_point(aes(y = Reef_mean, x = Year, color = "Reef level")) +
        geom_line(aes(y = y, x = Year, color = "Site level")) +
        geom_ribbon(aes(
          y = y, x = Year, ymin = ymin, ymax = ymax,
          fill = "Site level"
        ), alpha = 0.2) +
        scale_y_continuous("Cover", labels = function(x) x * 100) +
        scale_color_manual("",
          breaks = c("Reef level", "Site level"),
          values = c("red", "blue"), limits = c("Reef level", "Site level")
        ) +
        scale_fill_manual("",
          breaks = c("Reef level", "Site level"),
          values = c("red", "blue"), limits = c("Reef level", "Site level")
        ) +
        theme_bw() +
        ggtitle(.y)
    }
  ))
g$G |> patchwork::wrap_plots(nrow = 3)

## ----end

## Finer sampling design components ===============================

## ---- locs_obs_parameters
Number_of_transects_per_site <- 5
Depths <- 2
Number_of_frames_per_transect <- 100
Points_per_frame <- 5


## Note, the following are on the link scale
hcc_site_sigma <- 0.5        # variability in Sites within Locations
hcc_transect_sigma <- 0.2    # variability in Transects within Sites
hcc_sigma <- 0.1             # random noise

sc_site_sigma <- 0.05        # variability in Sites within Locations
sc_transect_sigma <- 0.02    # variability in Transects within Sites
sc_sigma <- 0.01             # random noise

ma_site_sigma <- 0.5        # variability in Sites within Locations
ma_transect_sigma <- 0.2    # variability in Transects within Sites
ma_sigma <- 0.1             # random noise

## ----end

## INLA ----------------------------------------------------

## Fixed locations ....................................................
## ---- fixed_locs_obs
set.seed(seed)
data_fixed_locs_obs <- data_fixed_locs_sf |>
  bind_cols(data_fixed_locs_sf |>
              st_coordinates() |>
              as.data.frame() |>
              dplyr::rename(Longitude = X, Latitude = Y)) |>
  st_drop_geometry() |>
  as.data.frame() |>
  group_by(Longitude, Latitude, Reef) |>
  crossing(
    Transect = paste0("T",1:Number_of_transects_per_site)) |>
  group_by(Site, .add = TRUE) |>
  mutate(
    SiteEffects_HCC = rnorm(1, 0, hcc_site_sigma),
    SiteEffects_SC = rnorm(1, 0, sc_site_sigma),
    SiteEffects_MA = rnorm(1, 0, ma_site_sigma)
  ) |>
  group_by(Transect, .add = TRUE) |>
  mutate(
    TransectEffects_HCC = rnorm(1, 0, hcc_transect_sigma),
    TransectEffects_SC = rnorm(1, 0, sc_transect_sigma),
    TransectEffects_MA = rnorm(1, 0, ma_transect_sigma)
  ) |>
  ungroup() |>
  mutate(
    HCC1 = HCC + SiteEffects_HCC +
      TransectEffects_HCC +
      rnorm(n(), 0, hcc_sigma),
    HCC2 = 100*plogis(HCC1),
    SC1 = SC + SiteEffects_SC + TransectEffects_SC +
      rnorm(n(), 0, sc_sigma),
    SC2 = 100*plogis(SC1),
    MA1 = MA + SiteEffects_MA + TransectEffects_MA
    + rnorm(n(), 0, ma_sigma),
    MA2 = 100*plogis(MA1)
  ) |>
  arrange(Reef, Site, Transect, Year) |>
  dplyr::select(Reef, Longitude, Latitude, Site,
    Transect, Year, HCC = HCC2, SC = SC2, MA = MA2) |>
  mutate(Year = 2021 - max(years) + Year,
    Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))
## ----end

## ---- fixed_locs_obs_depths
## The following are on a fold scale.
## Hence a value of 0.8, indicates that 
Depth_effect_multiplier <- 2

data_fixed_locs_obs <-
  data_fixed_locs_obs |>
  tidyr::crossing(Depth = seq(3, 10, length = Depths)) |>
  pivot_longer(cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value") |>
  group_by(Reef, Site, Transect, Year, Date) |>
  mutate(Value = Value + rev(sort(Depth_effect_multiplier *
                                    scale(rnorm(Depths))))) |>
  ungroup()
data_fixed_locs_obs |> head()
## ----end

## ---- fixed_locs_obs_fortify_data
## Need to split the percentage cover into point and frames
data_fixed_locs_obs <- data_fixed_locs_obs |>
  group_by(Reef,Site,Transect,Year,Depth,Date) |>
  mutate(Points = round(Number_of_frames_per_transect *
                          Points_per_frame *
                          (Value/sum(Value)),0),
    Points = ifelse(Points<0, 0, Points)) |>
  tidyr::uncount(Points) |>
  sample_n(n(), replace=FALSE) |>
  mutate(POINT_NO = rep_len(1:Points_per_frame, length = n()),
    ## FRAME = 1 + cumsum(POINT_NO) %/% (sum(1:Points_per_frame) + 1e-10)) |>
    FRAME = rep(1:Number_of_frames_per_transect, each=Points_per_frame, length = n())) |>
  ungroup() 

## a |> group_by(Reef, Site, Transect, Year, Depth, Group) |>
##     summarise(Count = n()) |>
##     ungroup(Group) |>
##     mutate(Total=sum(Count),
##            Cover = Count/Total)

reef_data_synthetic_fixed <- 
  data_fixed_locs_obs |> 
  mutate(PCODE = "SYNTHETIC-fixed",
    ID = 1:n(),
    CRUISE_CODE = paste0("SYNTHETIC",Year),
    REEF_NAME = Reef,
    AIMS_REEF_NAME = Reef,
    SECTOR = "synthetic",
    LATITUDE = Latitude,
    LONGITUDE = Longitude,
    SITE_NO = Site,
    TRANSECT_NO = Transect,
    SITE_DEPTH = Depth,
    REEF_ZONE = "-",
    REPORT_YEAR = Year,
    SURVEY_DATE = Date,
    FRAME = paste0(PCODE, "/", REEF_NAME, "/",
      REEF_ZONE, "/", SITE_NO, "/", SITE_DEPTH,
      "/", TRANSECT_NO, "/", REPORT_YEAR, "/", FRAME),
    POINT_NO = POINT_NO,
    FAMILY = NA,
    GROUP_DESC = Group,
    REEFPAGE_CATEGORY = paste0(Group,"_alt")
  ) |>
  dplyr::select(PCODE, ID, CRUISE_CODE, REEF_NAME,
    AIMS_REEF_NAME, SECTOR,
    LATITUDE, LONGITUDE, SITE_NO, TRANSECT_NO, SITE_DEPTH,
    REEF_ZONE, REPORT_YEAR, SURVEY_DATE, FRAME, POINT_NO,
    FAMILY, GROUP_DESC, REEFPAGE_CATEGORY)

write_csv(reef_data_synthetic_fixed,
  file = "../data/reef_data_synthetic_fixed.csv"
)
rmarkdown::paged_table(reef_data_synthetic_fixed |> head()) 

## ----end


## ---- hierarchical_schematic_fixed
#require(igraph)
#require(ggraph)
data <- reef_data_synthetic_fixed
 
edges <- data |>
  mutate(Origin = 'Root', Height = 6, Name = SECTOR, Level = 'Region',
    forground_colour = NA, background_colour = "white", point_size = "point_1") |>
  select(from = Origin, to = SECTOR, Height, Name, Level,
    forground_colour, background_colour, point_size) |>
  distinct() |>
  rbind(
    data |>
      mutate(Height = 5,
        Name = REEF_NAME,
        Level = 'Reef',
        forground_colour = NA,
        background_colour = "white",
        point_size = 'point_1') |>
      select(from = SECTOR,
        to = REEF_NAME,
        Height,
        Name,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(Height = 3,
        Name = SITE_NO,
        Level = 'Site',
        Site = paste(REEF_NAME, SITE_NO),
        forground_colour = NA,
        background_colour = "white",
        point_size = 'point_1') |>
      select(from = REEF_NAME,
        to = Site,
        Height,
        Name =  Site,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(
        Height = 2,
        Name = TRANSECT_NO,
        Level = "Transect",
        Site = paste(REEF_NAME, SITE_NO),
        Transect = paste(REEF_NAME, SITE_NO, TRANSECT_NO),
        forground_colour = NA,
        background_colour = "white",
        point_size = 'point_1') |> 
      select(from = Site,
        to = Transect,
        Height,
        Name = Transect,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(Height = 1,
        REPORT_YEAR = factor(REPORT_YEAR),
        Name = REPORT_YEAR,
        Level = 'Year',
        Site = paste(REEF_NAME, SITE_NO),
        Transect = paste(REEF_NAME, SITE_NO, TRANSECT_NO),
        Year = paste(REEF_NAME, SITE_NO, TRANSECT_NO, REPORT_YEAR),
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = "point_0") |>
      select(from = Transect,
        to = Year,
        Height,
        Name = Year,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(Height = 0,
        REPORT_YEAR = factor(REPORT_YEAR),
        ## Name = PHOTO,
        Level = 'Photo',
        Site = paste(REEF_NAME, SITE_NO),
        Transect = paste(REEF_NAME, SITE_NO, TRANSECT_NO),
        Year = paste(REEF_NAME, SITE_NO, TRANSECT_NO, REPORT_YEAR),
        Photo = paste(Year, str_extract(FRAME, "\\d+$")), 
        ## forground_colour = REPORT_YEAR,
        forground_colour = REPORT_YEAR,
        background_colour = NA,
        point_size = "point_0") |>
      select(from = Year,
        to = Photo,
        Height,
        Name = Photo,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  )  

reef_names <- data |>
  pull(REEF_NAME) |>
  unique()

edges <- edges |>
  filter(!(Level == "Photo" & !str_detect(from, paste(reef_names[1], "S1"))))

vertices <- edges |>
  select(Name, Height, Level, forground_colour, background_colour, point_size) |>
  distinct() |>
  rbind(data.frame(Name = 'Root', Height = 7, Level = 'Root',
    forground_colour = NA, background_colour = "white", point_size = 1)) |>
  mutate(
    Name2 = case_when(
      Level == 'Photo' ~ str_replace_all(Name, '.*',''),
      Level == 'Reef' ~ ifelse(Name == reef_names[1], Name, ""), #str_replace_all(Name, '.*',''),
      Level == 'Site' ~ ifelse(str_detect(Name, reef_names[1]), str_replace(Name, ".*(S.*)", "\\1"), ""),
      Level == 'Transect' ~ ifelse(str_detect(Name, paste(reef_names[1], "S1")), str_replace(Name, ".*(T.*)", "\\1"), ""),
      Level == 'Year' ~ ifelse(str_detect(Name, paste(reef_names[1], "S1")), str_replace(Name, ".*T.(.*)", "\\1"), ""),
      Level == 'Zone' ~ str_replace(Name, '[^_]*_[^_]*_', ''),
      !Level %in% c('Site') ~ str_replace(Name, "(.*)", "\\1")),
    size = Height
  )  

vertices <- vertices |>
  mutate(Name2 = ifelse(Level == "Colony", NA, Name2)) |>
  mutate(size = ifelse(str_detect(Name, ".*D$") & Level == "Zone", size - 0.2, size)) |>
  mutate(size = ifelse(str_detect(Name, ".*S$") & Level == "Zone", size + 0.2, size)) |>
  mutate(background_colour = ifelse(Name2 == "", NA, background_colour))


levels <- edges |> dplyr::select(Level, Height) |> distinct()
graph <- igraph::graph_from_data_frame(edges, vertices = vertices)

ggraph::ggraph(graph, layout = "dendrogram", circular = TRUE, height = size, repel = TRUE, ratio = 0.001) +
  ggraph::geom_edge_diagonal(aes(colour = edges$forground_colour), show.legend = FALSE) +
  ggraph::geom_node_point(aes(size = factor(point_size)), show.legend = FALSE) +
  ggraph::geom_node_label(aes(label = Name2, size = Level, 
    fill = background_colour),
    label.size = NA, repel = FALSE, show.legend = FALSE) +
  theme_void() +
  ## scale_y_continuous('', breaks = levels$Height, labels = levels$Level) +
  scale_size_manual(breaks = c("point_0", "point_1",'Photo','Reef', 'Region', 'Root', 'Site', 'Transect', 'Year', "A"),
    values = c(0, 1, 3,5,4,4,4,3,3, 0)) +
  scale_fill_manual(values = c("white"), na.value = "#00000000") 

detach("package:igraph")
detach("package:ggraph")
## ----end




## Random Locations ....................................................
## ---- random_locs_obs
set.seed(seed)
 
data_random_locs_obs <- data_random_locs_sf |>
  bind_cols(data_random_locs_sf |>
              st_coordinates() |>
              as.data.frame() |>
              dplyr::rename(Longitude = X, Latitude = Y)) |>
  st_drop_geometry() |>
  as.data.frame() |>
  group_by(Longitude, Latitude, Reef) |>
  crossing(Transect = paste0("T",
    1:Number_of_transects_per_site)) |>
  group_by(Site, .add = TRUE) |>
  mutate(
    SiteEffects_HCC = rnorm(1, 0, hcc_site_sigma),
    SiteEffects_SC = rnorm(1, 0, sc_site_sigma),
    SiteEffects_MA = rnorm(1, 0, ma_site_sigma)
  ) |>
  group_by(Transect, .add = TRUE) |>
  mutate(
    TransectEffects_HCC = rnorm(1, 0, hcc_transect_sigma),
    TransectEffects_SC = rnorm(1, 0, sc_transect_sigma),
    TransectEffects_MA = rnorm(1, 0, ma_transect_sigma)
  ) |>
  ungroup() |>
  mutate(
    HCC1 = HCC + SiteEffects_HCC + TransectEffects_HCC +
      rnorm(n(), 0, hcc_sigma),
    HCC2 = 100*plogis(HCC1),
    SC1 = SC + SiteEffects_SC + TransectEffects_SC +
      rnorm(n(), 0, sc_sigma),
    SC2 = 100*plogis(SC1),
    MA1 = MA + SiteEffects_MA + TransectEffects_MA +
      rnorm(n(), 0, ma_sigma),
    MA2 = 100*plogis(MA1)
  ) |>
  arrange(Reef, Site, Transect, Year) |>
  dplyr::select(Reef, Longitude, Latitude, Site, Transect,
    Year, HCC = HCC2, SC = SC2, MA = MA2) |>
  mutate(Year = 2021 - max(years) + Year,
    Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))
## ----end

## ---- random_locs_obs_depths
## The following are on a fold scale.
## Hence a value of 0.8, indicates that 
Depth_effect_multiplier <- 2

data_random_locs_obs <-
  data_random_locs_obs |>
  tidyr::crossing(Depth = seq(3, 10, length = Depths)) |>
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value"
  ) |>
  group_by(Reef, Site, Transect, Year, Date) |>
  mutate(Value = Value +
    rev(sort(Depth_effect_multiplier *
      scale(rnorm(Depths))))) |>
  ungroup()
## ----end

## ---- random_locs_obs_fortify_data
## Need to split the percentage cover into point and frames
data_random_locs_obs <- data_random_locs_obs |>
  group_by(Reef, Site, Transect, Year, Depth, Date) |>
  mutate(
    Points = round(Number_of_frames_per_transect *
      Points_per_frame * (Value / sum(Value)), 0),
    Points = ifelse(Points < 0, 0, Points)
  ) |>
  tidyr::uncount(Points) |>
  sample_n(n(), replace = FALSE) |>
  mutate(
    POINT_NO = rep_len(1:Points_per_frame, length = n()),
    ## FRAME = 1 + cumsum(POINT_NO) %/% (sum(1:Points_per_frame) + 1e-10)) |>
    FRAME = rep(1:Number_of_frames_per_transect,
      each = Points_per_frame, length = n()
    )
  ) |>
  ungroup()

reef_data_synthetic_random <-
  data_random_locs_obs |>
  mutate(
    PCODE = "SYNTHETIC-random",
    ID = 1:n(),
    CRUISE_CODE = paste0("SYNTHETIC", Year),
    REEF_NAME = Reef,
    AIMS_REEF_NAME = Reef,
    SECTOR = "synthetic",
    LATITUDE = Latitude,
    LONGITUDE = Longitude,
    SITE_NO = Site,
    TRANSECT_NO = Transect,
    SITE_DEPTH = Depth,
    REEF_ZONE = "-",
    REPORT_YEAR = Year,
    SURVEY_DATE = Date,
    FRAME = paste0(PCODE, "/", REEF_NAME, "/", REEF_ZONE,
      "/", SITE_NO, "/", SITE_DEPTH, "/", TRANSECT_NO,
      "/", REPORT_YEAR, "/", FRAME),
    POINT_NO = POINT_NO,
    FAMILY = NA,
    GROUP_DESC = Group,
    REEFPAGE_CATEGORY = paste0(Group, "_alt")
  ) |>
  dplyr::select(
    PCODE, ID, CRUISE_CODE, REEF_NAME, AIMS_REEF_NAME, SECTOR,
    LATITUDE, LONGITUDE, SITE_NO, TRANSECT_NO, SITE_DEPTH,
    REEF_ZONE, REPORT_YEAR, SURVEY_DATE, FRAME, POINT_NO,
    FAMILY, GROUP_DESC, REEFPAGE_CATEGORY
  )

write_csv(reef_data_synthetic_random, file = "../data/reef_data_synthetic_random.csv")
rmarkdown::paged_table(reef_data_synthetic_random |> head())
## ----end


## ---- random_download_data_btn
reef_data_synthetic_random |>
  download_this(
    output_name = "reef data synthetic random",
    output_extension = ".csv",
    button_label = "Download data",
    button_type = "warning",
    has_icon = TRUE,
    icon = "fa fa-save"
  )
## ----end


## ---- hierarchical_schematic_random
#require(igraph)
#require(ggraph)
data <- reef_data_synthetic_random

first_reef <- data |>
  arrange(REPORT_YEAR, REEF_NAME) |> 
  dplyr::filter(REPORT_YEAR == first(REPORT_YEAR)) |>
  droplevels() |>
  dplyr::filter(REEF_NAME == REEF_NAME[1]) |>
  pull(REEF_NAME) |>
  unique()
data <- data |>
  nest_by(REPORT_YEAR, REEF_NAME, .keep = TRUE) |>
  mutate(data = list({
    data |>
      filter(!(REPORT_YEAR != 2010 & SITE_NO == last(SITE_NO))) |>
      mutate(SITE_NO = ifelse(REPORT_YEAR != 2010, NA, SITE_NO)) |>
      mutate(Site = paste(REEF_NAME, SITE_NO)) |> 
      mutate(TRANSECT_NO = ifelse(REPORT_YEAR == 2010 &
        REEF_NAME == first_reef,
        #Site == first(Site),
        TRANSECT_NO, NA
      )) |> 
      mutate(FRAME = ifelse(REPORT_YEAR == 2010 &
                              REEF_NAME == first_reef &
       SITE_NO == "S1",# & TRANSECT_NO == "T1",
        #Site == first(Site),
        FRAME, NA
      ))
  })) |>
    ungroup() |> 
  dplyr::select(-REPORT_YEAR, -REEF_NAME) |> 
  unnest(data) 

edges <- data |>
  mutate(Origin = 'Root', Height = 6, Name = SECTOR, Level = 'Region',
    forground_colour = NA, background_colour = "white", point_size = "point_1") |>
  select(from = Origin, to = SECTOR, Height, Name, Level,
    forground_colour, background_colour, point_size) |>
  distinct() |>
  rbind(
    data |>
      mutate(Height = 5,
        REPORT_YEAR = factor(REPORT_YEAR),
        Name = REPORT_YEAR,
        Level = 'Year',
        Year = paste(REPORT_YEAR),
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = "point_0") |>
      select(from = SECTOR,
        to = Year,
        Height,
        Name = Year,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(Height = 4,
        Name = REEF_NAME,
        Level = 'Reef',
        Year = paste(REPORT_YEAR),
        Reef = paste(Year, REEF_NAME),
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = 'point_1') |>
      select(from = Year,
        to = Reef,
        Height,
        Name =  Reef,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct() 
       )  |> 
  rbind(
    data |>
      mutate(Height = 3,
        Name = SITE_NO,
        Level = 'Site',
        Year = paste(REPORT_YEAR),
        Reef = paste(Year, REEF_NAME),
        Site = ifelse(is.na(SITE_NO), NA, paste(Reef, SITE_NO)),
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = 'point_1') |>
      select(from = Reef,
        to = Site,
        Height,
        Name =  Site,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(
        Height = 2,
        Name = TRANSECT_NO,
        Level = "Transect",
        Year = paste(REPORT_YEAR),
        Reef = paste(Year, REEF_NAME),
        Site = ifelse(is.na(SITE_NO), NA, paste(Reef, SITE_NO)),
        Transect = ifelse(is.na(TRANSECT_NO), NA, paste(Site, TRANSECT_NO)),
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = 'point_1') |> 
      select(from = Site,
        to = Transect,
        Height,
        Name = Transect,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  ) |> 
  rbind(
    data |>
      mutate(Height = 0,
        REPORT_YEAR = factor(REPORT_YEAR),
        ## Name = PHOTO,
        Level = 'Photo',
        Year = paste(REPORT_YEAR),
        Reef = paste(Year, REEF_NAME),
        Site = ifelse(is.na(SITE_NO), NA, paste(Reef, SITE_NO)),
        Transect = ifelse(is.na(TRANSECT_NO), NA, paste(Site, TRANSECT_NO)),
        Photo = paste(Transect, str_extract(FRAME, "\\d+$")), 
        ## forground_colour = REPORT_YEAR,
        forground_colour = REPORT_YEAR,
        background_colour = "white",
        point_size = "point_0") |>
      select(from = Transect,
        to = Photo,
        Height,
        Name = Photo,
        Level,
        forground_colour,
        background_colour,
        point_size) |>
      distinct()
  )  

reef_names <- data |>
  pull(REEF_NAME) |>
  unique()

edges <- edges |>
  filter(!is.na(to)) |>
  filter(!str_detect(to, "NA"))
##   ##filter(!(Level == "Photo" & !str_detect(from, paste("2011")))) |> 
##   filter(!(Level == "Site" & !str_detect(from, paste("2011"))))

edges <- edges |>
  arrange(from, to)

vertices <- edges |>
  select(Name, Height, Level, forground_colour, background_colour, point_size) |>
  distinct() |>
  rbind(data.frame(Name = 'Root', Height = 7, Level = 'Root',
    forground_colour = NA, background_colour = "white", point_size = 1)) |>
  mutate(
    Name2 = case_when(
      ## Level == 'Photo' ~ str_replace_all(Name, '.*',''),
      # Level == 'Reef' ~ str_replace_all(Name, '[0-9]{4} (.*)','\\1'),
      ## Level == 'Reef' ~ ifelse(Name == reef_names[1], Name, ""), #str_replace_all(Name, '.*',''),
      Level == "Reef" ~ ifelse(str_detect(Name, paste("2010")) &
        str_detect(Name, first_reef),
      str_replace(Name, "[0-9]{4} (.*)", "\\1"), ""
      ),
      Level == "Site" ~ ifelse(str_detect(Name, paste("2010")) &
        str_detect(Name, first_reef),
      str_replace(Name, ".*(S.*)", "\\1"), ""
      ),
      Level == "Transect" ~ ifelse(str_detect(Name, paste("2010")) &
        str_detect(Name, "S1"),
      str_replace(Name, ".*(T.*)", "\\1"), ""
      ),
      ## Level == "Photo" ~ ifelse(str_detect(Name, paste("2010")) &
      ##                             str_detect(Name, "S1") &
      ##                             str_detect(Name, "T1"),
      ## str_replace(Name, ".* (\\d+$)", "\\1"), ""
      ## ),
      Level == "Photo" ~ "",
      !Level %in% c("Site") ~ str_replace(Name, "(.*)", "\\1")
    ),
    size = Height
  ) |>
  filter(!str_detect(Name, "NA"))

vertices <- vertices |>
  mutate(Name2 = ifelse(Level == "Colony", NA, Name2)) |>
  mutate(size = ifelse(str_detect(Name, ".*D$") & Level == "Zone", size - 0.2, size)) |>
  mutate(size = ifelse(str_detect(Name, ".*S$") & Level == "Zone", size + 0.2, size)) |>
  mutate(background_colour = ifelse(Name2 == "", NA, background_colour))


levels <- edges |> dplyr::select(Level, Height) |> distinct()
graph <- igraph::graph_from_data_frame(edges, vertices = vertices)

ggraph::ggraph(graph, layout = "dendrogram", circular = TRUE, height = size, repel = TRUE, ratio = 1000) +
  ggraph::geom_edge_diagonal(aes(colour = forground_colour), show.legend = FALSE) +
  ggraph::geom_node_point(aes(size = factor(point_size)), show.legend = FALSE) +
  ggraph::geom_node_label(aes(label = Name2, size = Level, 
    fill = background_colour),
    label.size = NA, repel = FALSE, show.legend = FALSE) +
  theme_void() +
    ggraph::scale_edge_colour_discrete() +
  ## scale_y_continuous('', breaks = levels$Height, labels = levels$Level) +
  scale_size_manual(breaks = c("point_0", "point_1",'Photo','Reef', 'Region', 'Root', 'Site', 'Transect', 'Year', "A"),
    values = c(0, 1, 3,3,4,4,4,3,3, 0)) +
  scale_fill_manual(values = c("white"), na.value = "#00000000") 

detach("package:igraph")
detach("package:ggraph")
## ----end
