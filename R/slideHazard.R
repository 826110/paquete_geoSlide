#' Calcular pendiente de un DEM
#'
#' Calcula la pendiente de un Modelo Digital de Elevacion(DEM) en grados.
#'
#'
#' @param dem Un objeto SpatRaster que represente el Modelo Digital de Elevacion
#' (DEM).
#'
#' @return Un objeto SpatRaster con valores de pendientes expresados en unidades
#' de grados.
#'
#' @export
#' @examples
#' library(terra)
#'
#' dem <- rast(nrows = 50, ncols = 50,
#'             xmin = 0, xmax = 1000,
#'             ymin = 0, ymax = 1000,
#'             crs = "EPSG:25830")
#'
#' values(dem) <- runif(ncell(dem), 0, 100)
#'
#' slope <- slope_func(dem)
#' plot(slope)
#'

slope_func <- function(dem) {

  #Check objeto de entrada#

  if (!inherits(dem, "SpatRaster")) {
    stop("dem ha de ser un objeto SpatRaster.")
  }

  #Check sistema de coordenadas (CRS)#

  if (is.na(terra::crs(dem))) {
    stop("dem debe de tener un sistema de cordenadas (CRS) definido")
  }

  #Calculo de pendientes#

  slope <- terra::terrain(dem,
                          v = "slope",
                          unit = "degrees")

  return(slope)
}

#' Distancia a fallas geologicas
#'
#' Calcula la distancia Euclidiana de cada celda del raster a la linea de falla
#' mas proxima
#'
#' @param dem Un objeto SpatRaster que sirve como referencia espacial
#' @param faults Un objeto sf con las diferentes fallas poligonizadas
#'
#' @return Un objeto SpatRaster con valores de distancia
#' @export
#'
#' @examples
#' library(terra)
#' library(sf)
#'
#' dem <- rast(nrows = 50, ncols = 50,
#'             xmin = 0, xmax = 1000,
#'             ymin = 0, ymax = 1000,
#'             crs = "EPSG:25830")
#'
#' values(dem) <- runif(ncell(dem), 0, 100)
#'
#' line <- st_linestring(matrix(c(200,200,
#'                                800,800),
#'                              ncol = 2,
#'                              byrow = TRUE))
#'
#' faults <- st_sf(
#'   geometry = st_sfc(line, crs = 25830)
#' )
#'
#' dist <- distance_to_fault(dem, faults)
#' plot(dist)

distance_to_fault <- function(dem, faults) {

 #Check de objetos de entrada y Sistema de Referencia#

  if (!inherits(dem, "SpatRaster")) {
    stop("dem debe de ser un objeto SpatRaster.")
  }

  if (!inherits(faults, "sf")) {
    stop("faults debe de ser un objeto sf.")
  }

  if (is.na(terra::crs(dem)) || is.na(sf::st_crs(faults))) {
    stop("Ambos objetos (dem y faults) deben de tener un CRS definido")
  }

  if (sf::st_crs(faults)$wkt != terra::crs(dem)) {
    stop("Ambos objetos (dem y faults) deben de tener el mismo CRS.")
  }

  # Convertir sf → SpatVector
  faults_vect <- terra::vect(faults)

  #Calcular distancia#

  dist_raster <- terra::distance(dem, faults_vect)

  return(dist_raster)
}

#' Normalizar los raster
#'
#' Normaliza un SpatRaster usando el método min-max.
#'
#' @param x Objetos SpatRaster heredados de las funciones distance_to_fault y
#' slope_func
#'
#' @return SpatRaster normalizado entre 0 y 1
#' @export
#' @examples
#' library(terra)
#'
#' r <- rast(nrows = 10, ncols = 10,
#'           xmin = 0, xmax = 100,
#'           ymin = 0, ymax = 100,
#'           crs = "EPSG:25830")
#'
#' values(r) <- runif(ncell(r), 10, 50)
#'
#' r_norm <- normalize_minmax(r)
#' terra::global(r_norm, c("min", "max"))

normalize_minmax <- function(x) {

  min_val <- terra::global(x, "min", na.rm = TRUE)[1,1]
  max_val <- terra::global(x, "max", na.rm = TRUE)[1,1]

  if (min_val == max_val) {
    stop("No se puede normalizar: valores constantes.")
  }

  x_norm <- (x - min_val) / (max_val - min_val)

  return(x_norm)
}

#' Indice de peligrosidad de deslizamiento
#'
#' Calcula un índice simple de peligrosidad por deslizamiento
#' combinando pendiente y distancia a fallas.
#'
#' @param slope SpatRaster de pendiente (normalizado 0–1) heredado de la funcion
#' slope_func.
#' @param distance SpatRaster de distancia a fallas (normalizado 0–1) heredado
#' de la funcion distance_to_fault.
#' @param w_slope Peso de la pendiente (por defecto = 0.5).
#' @param w_distance Peso de la distancia (por defecto = 0.5).
#'
#' @return SpatRaster con índice de peligrosidad.
#' @export
#' @examples
#' library(terra)
#'
#' # Crear raster simulado
#' r <- rast(nrows = 50, ncols = 50,
#'           xmin = 0, xmax = 1000,
#'           ymin = 0, ymax = 1000,
#'           crs = "EPSG:25830")
#'
#' values(r) <- runif(ncell(r), 0, 100)
#'
#' slope_n <- normalize_minmax(r)
#' dist_n  <- normalize_minmax(r)
#'
#' hazard <- hazard_index(slope_n, dist_n,
#'                             w_slope = 0.6,
#'                             w_distance = 0.4)
#'
#' plot(hazard)

hazard_index <- function(
    slope,
    distance,
    w_slope = 0.5,
    w_distance = 0.5)
  {

  # Verificar pesos
  if ((w_slope + w_distance) != 1) {
    stop("Los pesos deben sumar 1.")
  }

  if (any(c(w_slope, w_distance) < 0)) {
    stop("Los pesos no pueden ser negativos.")
  }

  # Invertir distancia (más cerca = más peligro)
  distance_inv <- 1 - distance

  hazard <- (slope * w_slope) +
    (distance_inv * w_distance)

  return(hazard)
}
