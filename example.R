source("unknown_pleasuR.R")
library(raster)
library(tigris)
library(tidycensus)
library(stringr)
library(purrr)
library(gstat)
library(dplyr)
library(tidyr)

vars <- list(
  # Show variables with...
  # vars <- load_variables(2019, "acs5", cache = TRUE)
  tot_rental_hh = "B25003_003", 
  black_renters = "B25003B_003", 
  latinx_renters = "B25003I_003"
)

var_list <- unname(unlist(vars))
var_names <- names(vars)
var_tib <- tibble(old=var_list, new=var_names)
table <- get_acs(
  geo = "tract", 
  state = "MA",
  year = 2020, 
  variables = var_list, 
  geometry = FALSE
)

tracts_geom <- tracts("MA", year = 2020, cb = TRUE)

tracts_table <- get_acs(
    geography="tract", 
    variables = var_list, 
    state="MA",
    year=2020, 
    survey='acs5',
    geometry=FALSE
  ) %>%
  rename_with(str_to_lower) %>%
  dplyr::select(-c('moe')) %>%
  pivot_wider(names_from = 'variable', values_from = 'estimate') %>%
  rename(set_names(var_tib$old, var_tib$new)) %>%
  rename(
    id = geoid
  ) %>%
  mutate(
    rent_pct_black = black_renters / tot_rental_hh, 
    rent_pct_latinx = latinx_renters / tot_rental_hh
  ) 

tracts <- tracts_geom %>%
  st_transform(2249) %>%
  st_make_valid() %>%
  rename_with(str_to_lower) %>%
  rename(
    id = geoid
  ) %>%
  left_join(tracts_table, by = "id")

dims <- get_dims(tracts, n = 100, type = "vertical")

raster_black_renters <- interpolate(
  raster(tracts, res=500),
  gstat(
    formula = rent_pct_black ~ 1, 
    nmax = 20, 
    set = list(idp = 2), 
    data = st_centroid(drop_na(tracts, rent_pct_black))
  )
) %>%
  mask(tracts)

raster_latinx_renters <- interpolate(
  raster(tracts, res=500),
  gstat(
    formula = rent_pct_latinx ~ 1, 
    nmax = 20, 
    set = list(idp = 2), 
    data = st_centroid(drop_na(tracts, rent_pct_latinx))
  )
) %>%
  mask(tracts)

lines <- tracts %>%
  st_union() %>%
  st_regular_lines(
    dims = dims,
    mask = TRUE
  )

black_renters <- lines %>%
  st_unknown_pleasures(
    raster_black_renters,
    dims = dims,
    sample_size = 250, 
    bleed_factor = 5,
    mode = "xyz",
    polygon = FALSE
  )



latinx_renters <- lines %>%
  st_unknown_pleasures(
    raster_latinx_renters, 
    dims = dims,
    sample_size = 250, 
    bleed_factor = 2,
    polygon = TRUE
  )