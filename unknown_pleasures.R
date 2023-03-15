get_axes <- function(type = "horizontal") {
  #' Get information identifying orientation of lines.
  #'
  #' @param type String indicating direction of lines---options are "horizontal" and "vertical"
  #' @returns A list of axes.
  #' @export
  if (type == "horizontal") {
    list("int_dim" = "y", "edge_dim" = "x")
  } else if (type == "vertical") {
    list("int_dim" = "x", "edge_dim" = "y")
  } else {
    stop("Invalid line type - expects horizontal or vertical.")
  }
}

get_interval <- function(bbox, n, type = "horizontal") {
  #' Generate list of interval-relevant values.
  #'
  #' @param bbox Object of type `bbox`.
  #' @param n Number of regularly-spaced lines to generate.
  #' @param type String indicating direction of lines---options are "horizontal" and "vertical".
  #' @returns A list of interval-relevant values.
  #' @export
  axes <- get_axes(type)
  int_min <- unname(bbox[paste0(axes$int_dim, "min")])
  int_max <- unname(bbox[paste0(axes$int_dim, "max")])
  list(
    "int_min" = int_min, 
    "int_max" = int_max, 
    "interval" = (int_max - int_min) / (n - 1)
    )
}

st_regular_lines <- function(df, n, mask = TRUE, type = "horizontal") {
  #' Generate n regularly-spaced lines over `df` extent.
  #'
  #' @param df An `sf` dataframe.
  #' @param n Number of regularly-spaced lines to generate.
  #' @param mask Whether resulting lines should be clipped to df extent.
  #' @param type String indicating direction of lines---options are "horizontal" and "vertical".
  #' @returns A list of interval-relevant values.
  #' @export
  print("Generating regularly spaced lines over input bbox...")
  bbox <- df %>%
    sf::st_bbox()
  axes <- get_axes(type)
  intervals <- get_interval(bbox, n = n)
  edge_min <- unname(bbox[paste0(axes$edge_dim, "min")])
  edge_max <- unname(bbox[paste0(axes$edge_dim, "max")])
  line_positions <- seq(
      from = intervals$int_min, 
      to = intervals$int_max, 
      by = intervals$interval
    )
  lines <- sf::st_sfc(crs = sf::st_crs(df))
  for (line in line_positions) {
    line_const <- sf::st_sfc(
      sf::st_linestring(
        matrix(
          c(edge_min, line, edge_max, line), 
          ncol = 2, 
          byrow = TRUE
          )
      ),
      crs = sf::st_crs(df)
    )
    lines <- sf::st_sfc(
        c(lines, line_const),
        crs = sf::st_crs(df)
    )
  }
  sf <- sf::st_as_sf(data.frame(id = 1:n, geometry = lines)) 
  if (mask) {
    print("Clipping lines to input layer extent...")
    sf %>%
      sf::st_intersection(
        sf::st_geometry(df)
      ) %>%
      dplyr::rowwise() %>%
      dplyr::filter(st_geometry_type(geometry) != "POINT") %>%
      dplyr::ungroup() %>% 
      sf::st_cast("MULTILINESTRING") %>%
      sf::st_cast("LINESTRING") %>%
      dplyr::mutate(
        id = dplyr::row_number()
      )
  } else {
    sf
  }
}

st_unknown_pleasures <- function(
    lines, 
    raster, 
    n,
    sample_size = 250,
    bleed_factor = 1.5,
    polygon = TRUE) {
  #' Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.
  #'
  #' @param lines A spatial dataframe containing regularly spaced lines.
  #' @param n Number of regularly-spaced lines to generate.
  #' @param sample_size Interval, in meters, to sample along lines.
  #' @param bleed_factor How much should maxima bleed into the area of the line above or below.
  #' @param polygon If `TRUE`, outputs polygons. If `FALSE`, outputs lines.
  #' @returns A spatial dataframe containing "Unknown Pleasures" features.
  #' @export
  print("Sampling along lines...")
  elevated_lines <- lines %>%
    dplyr::mutate(
      geometry = sf::st_line_sample(
        geometry, 
        density = (1 / units::as_units(sample_size, "m")),
        type = "regular"
      )
    ) %>%
    sf::st_cast("POINT") 
  print("Extracting elevations at sample points...")
  elevated_lines$elev <- raster::extract(raster, elevated_lines)
  max <- max(abs(elevated_lines$elev), na.rm = TRUE)
  ints <- get_interval(sf::st_bbox(elevated_lines), n = n)
  scale <- ints$interval / max
  print("Performing Affine Transform on points...")
  elevated_lines <- elevated_lines %>%
    tidyr::drop_na(elev) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      elev_scaled = elev * (scale * bleed_factor)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      geometry = geometry + c(0, elev_scaled)
    ) %>%
    dplyr::ungroup() %>%
    sf::st_set_crs(sf::st_crs(lines)) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize() %>%
    sf::st_cast("LINESTRING") %>%
    dplyr::ungroup()
  
  if (polygon) {
    print("Building polygons...")
    lines <- lines %>%
      dplyr::filter(id %in% dplyr::pull(elevated_lines, id))
    dplyr::bind_rows(lines, elevated_lines) %>%
      sf::st_cast("POINT") %>%
      dplyr::group_by(id) %>%
      dplyr::summarize() %>%
      sf::st_cast("POLYGON") %>%
      dplyr::arrange(desc(id)) %>%
      dplyr::ungroup ()
  }
}