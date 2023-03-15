source("unknown_pleasures.R")
library(raster)
library(tidycensus)
library(stringr)
library(purrr)
library(gstat)

vars <- list(
  # Show variables with...
  # vars <- load_variables(2019, "acs5", cache = TRUE)
  tot_hh = "B25003_001",
  tot_rental_hh = "B25003_003", 
  white_renters = "B25003H_003", 
  black_renters = "B25003B_003", 
  indig_renters = "B25003C_003",
  asian_renters = "B25003D_003", 
  pi_renters = "B25003E_003",
  other_renters = "B25003F_003",
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

library(tigris)
tracts_geom <- tracts("MA", year = 2020, cb = TRUE)

tracts_table <- get_acs(geography="tract", 
                        variables = var_list, 
                        state="MA",
                        year=2020, 
                        survey='acs5',
                        geometry=FALSE) %>%
  rename_with(str_to_lower) %>%
  dplyr::select(-c('moe')) %>%
  pivot_wider(names_from = 'variable', values_from = 'estimate') %>%
  rename(set_names(var_tib$old, var_tib$new)) %>%
  rename(
    id = geoid
  ) %>%
  mutate(
    rent_pct = tot_rental_hh / tot_hh,
    rent_pct_white = white_renters / tot_rental_hh,
    rent_pct_nonwhite = 1 - rent_pct_white,
    rent_pct_black = black_renters / tot_rental_hh, 
    rent_pct_asian = asian_renters / tot_rental_hh,
    rent_pct_aapi = (asian_renters + pi_renters) / tot_rental_hh,
    rent_pct_latinx = latinx_renters / tot_rental_hh,
    rent_pct_indig = indig_renters / tot_rental_hh,
    rent_pct_other = other_renters / tot_rental_hh
  ) 

tracts <- tracts_geom %>%
  st_transform(2249) %>%
  st_make_valid() %>%
  rename_with(str_to_lower) %>%
  rename(
    id = geoid
  ) %>%
  left_join(tracts_table, by = "id")

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
    n = 100,
    mask = TRUE,
    type = "horizontal"
  )

black_renters <- lines %>%
  st_unknown_pleasures(
    raster_black_renters, 
    n = 100, 
    sample_size = 500, 
    bleed_factor = 2
  )

latinx_renters <- lines %>%
  st_unknown_pleasures(
    raster_latinx_renters, 
    n = 100, 
    sample_size = 500, 
    bleed_factor = 2
  )