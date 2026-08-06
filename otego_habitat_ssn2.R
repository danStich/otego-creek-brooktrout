# Load Libraries and Define Cores ----
# . Libraries ----
library(SSN2)
library(SSNbler)
library(sf)
library(tidyverse)
library(parallel)

# . Parallel settings ----
detectCores() # How many cores do I have?
num.cores <- 8  # Don't use all the cores. Leave a few to do other things.

# Read data ----
# . Paths for lsn and ssn files ----
ssn_path1<- "GIS/HSI"
lsn_path3 <- "GIS/LSN3" # Pre-processing flowline result files
lsn_path4 <- "GIS/LSN4" # Empty dir for processing that follows

# . Edges ----
edges3 <- st_read(paste0(lsn_path3, "/otego_creek.shp"))

# .. Re-check topology for remaining errors ----
edges <- lines_to_lsn(
  streams = edges3,
  lsn_path = lsn_path4,
  snap_tolerance = 1,
  check_topology = TRUE,
  topo_tolerance = 20,
  overwrite = TRUE,
  verbose = TRUE,
  remove_ZM = TRUE,
  use_parallel = TRUE,
  no_cores = num.cores)


# . Observed and predicted points ----
# .. Observed sites ----
obs_points <- st_read("GIS/sites/bkt-hsi-elevations.shp")

# ... Renaming columns ----
names(obs_points) <- gsub(pattern = "\\.", 
                          replacement = "",
                          x = names(obs_points))


# ... Calculating HSIs ----
obs_points <- obs_points %>%
  mutate(across(
    .cols = 20:33,  # only numeric columns
    ~ if_else(
      is.na(.),
      mean(., na.rm = TRUE),  # column mean
      .
    )
  ))
                     
# Variables
v1 <- obs_points$V1
v3 <- obs_points$V3
v4 <- obs_points$V4
v6A <- obs_points$V6A
v6B <- obs_points$V6B
v8 <- obs_points$V8
v9 <- obs_points$V9
v10 <- obs_points$V10
v11 <- obs_points$V11
v12 <- obs_points$V12
v13 <- obs_points$V13
v14 <- 1
v15 <- obs_points$V15
v16 <- obs_points$V16
v17 <- obs_points$V17


# .... Other ----
CO <- ((((v9*v16)^(1/2) + v11)/2) * 
         (v1*v3*v12*v13*v14*v17)^(1/6))^(1/2)

# .... Adult ----
CA <- vector(mode = "numeric", length = length(v4))

CA[v6A > (v10 * v15)^(1/2)] <- ((v4 * v6A * ((v10 *v15)^(1/2))^(1/3)))[v6A > (v10 * v15)^(1/2)]
CA[v6A < (v10 * v15)^(1/2)] <- ((v4 * ((v10 *v15)^(1/2))^(1/2)))[v6A < (v10 * v15)^(1/2)]
CA[(v4 <= 0.4) | ((v10 * v15)^(1/2) <= 0.4) ] <- min(CA)
CA <- (CA * CO) ^ (1/2)

# .... Juvenile ----
CJ <- (v6B + v10 + v15)/3
CJ[(v6B <= 0.4) | (v10 <= 0.4) | (v15 <= 0.4) ] <- min(CJ)
CJ <- (CJ * CO) ^ (1/2)

# .... Fry ----
CF <- (v10*((v8 * v16)^(1/2)))^(1/2)
CF[(v10 <= 0.4) | ((v8 * v16)^(1/2) <= 0.4) ] <- min(CF)
CF <- (CF * CO) ^ (1/2)

# .... Combined HSI ----
# Formula
HSI <- (CA * CJ * CF * CO)^(1/4)

# Conditions
for(i in 1:length(HSI)){
  if(CA[i] <= 0.4 | CJ[i] <= 0.4 | CF[i] <= 0.4 ){
    HSI[i] <- min(CA[i], CJ[i], CF[i])
  }
}

HSI[CA <= HSI] <- CA


# .... Add HSIs to the obs_points ----
obs_points$CA <- CA
obs_points$CJ <- CJ
obs_points$CF <- CF
obs_points$HSI <- HSI


# .. Prediction sites ----
preds <- st_read("GIS/mid_atlantic_preds/otego_preds_muahaha.shp")


# Prepare LSN (landscape network) ----
# . Add observation sites to lsn ----
obs <- sites_to_lsn(
  sites = obs_points,
  edges = edges,
  snap_tolerance = 200,
  save_local = TRUE,
  lsn_path = lsn_path4,
  file_name = "sites.gpkg",
  overwrite = TRUE
)


# . Add prediction points to lsn ----
bkt_preds <- sites_to_lsn(
  sites = preds,
  edges = edges,
  snap_tolerance = 200,
  save_local = TRUE,
  lsn_path = lsn_path4,
  file_name = "preds.gpkg",
  overwrite = TRUE
)


# . Generate upstream distance for edges ----
edges <- updist_edges(
  edges = edges,
  lsn_path = lsn_path4,
  calc_length = TRUE,
  length_col = "Length",
  save_local = TRUE,
  overwrite = TRUE,
)


# . Generate the sites list ---- 
# This is a named list
sitelist = list(obs = obs, bkt_preds = bkt_preds)


# . Compute updist on sites and preds ----
site.list <- updist_sites(
  sites = sitelist,
  edges = edges,
  length_col = "Length",
  lsn_path = lsn_path4,
  save_local = TRUE,
  overwrite = TRUE
)


# Generate additive function values ----
# . Edges ----
edges <- afv_edges(
  edges = edges,
  lsn_path = lsn_path4,
  infl_col = "TotDASqKM",
  # infl_col = "CUMDRAINAG",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  save_local = TRUE,
  overwrite = TRUE
)


# . Sites ----
# using the site list created above
site.list <- afv_sites(
  sites = site.list,
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = lsn_path4,
  overwrite = TRUE
)

# Do these need to be given a non-zero afv or fixed??
# edges$afvArea[edges$afvArea == 0] <- .01

# Generate the SSN object ----
# Assemble ssn
covariate_ssn <- ssn_assemble(
  edges = edges,
  lsn_path = lsn_path4,
  obs_sites = site.list$obs, # This one is a df
  preds_list = site.list["bkt_preds"], # This one is a list
  ssn_path = ssn_path1,
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)

# Generate hydrologic distance matrices ----
# Calculate distances between features
ssn_create_distmat(
  ssn.object = covariate_ssn,
  predpts = c("bkt_preds"),
  among_predpts = TRUE,
  overwrite = TRUE
)

# Fit a test model ----
# SAMPLE_1 is elevation
# Method is eDNA or Electrofishinghttp://127.0.0.1:11869/graphics/49df2473-481b-4516-8e41-22039ac9ef15.png
ssn_mod <- ssn_glm(
  formula = CA ~ Elevation1,
  ssn.object = covariate_ssn,
  family = "gaussian",
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)

# . Summary ----
summary(ssn_mod)


# . Plot the residuals ----
# plot(ssn_mod, which = 1)


# . Predictions ----
aug_preds <- augment(ssn_mod, newdata = "bkt_preds",
                     type.predict = "response")

# ... Plot ----
plotter <- aug_preds %>% 
  arrange(.fitted)

bkt_plot <- ggplot() +
  geom_sf(data = plotter, aes(color = .fitted), 
          size = 3, alpha = 0.5) +
  geom_sf(data = covariate_ssn$edges) +
  geom_sf(data = obs_points) +
  scale_color_viridis_c(limits = c(0, 1), option = "H") +
  labs(color = "") +
  ylab("Latitude") +
  xlab("Longitude") +
  theme_bw()

bkt_plot

# Write plot to file 
jpeg(filename = "CA_preds.jpg",
     res = 300,
     height = 2400,
     width = 1800)
bkt_plot
dev.off()
