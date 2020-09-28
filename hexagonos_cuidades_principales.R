rm(list=ls())
gc()

######################Lectura de Librerías###############
library(dplyr)
library(tidyr)
library(sp)
library(raster)
library(rgeos)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
library(geojsonsf)
library(rgdal)
library(geojson)
library(mapview)
library(geojsonio)
library(sf)
########################################################3

# leemos el shape politico de Colombia
shp_pol <- st_read("otras cuidades/POLITICO/MGN_MPIO_POLITICO.shp")
#mapview(shp_pol)

# Seleccionamos los shapes de las ciudades principales
cali <- subset(shp_pol,DPTO_CCDGO=="76" & MPIO_CCDGO=="001") 
medellin <- subset(shp_pol,DPTO_CCDGO=="05" & MPIO_CCDGO=="001")
cartagena <- subset(shp_pol,DPTO_CCDGO=="13" & MPIO_CCDGO=="001")
barranquilla <- subset(shp_pol,DPTO_CCDGO=="08" & MPIO_CCDGO=="001")
bucaramanga <- subset(shp_pol,DPTO_CCDGO=="68" & MPIO_CCDGO=="001")

#############################################################

# los convertimos en geo_json y los guardamos
cali <- as.geojson(cali)
write(cali, "otras cuidades/ciudades/cali.geojson")
medellin <- as.geojson(medellin)
write(medellin, "otras cuidades/ciudades/medellin.geojson")
cartagena <- as.geojson(cartagena)
write(cartagena, "otras cuidades/ciudades/cartagena.geojson")
barranquilla <- as.geojson(barranquilla)
write(barranquilla, "otras cuidades/ciudades/barranquilla.geojson")
bucaramanga <- as.geojson(bucaramanga)
write(bucaramanga, "otras cuidades/ciudades/bucaramanga.geojson")

#############################################################

# leemos la funcion para los hexagonos

# Función para subdivir un shape en hexagonos#

make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    x <- gBuffer(x, byid=TRUE, width=0)
    g <- gBuffer(g, byid=TRUE, width=0)
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

###########################################################

# leemos los geojson 
cali <- geojson_read( "otras cuidades/ciudades/cali.geojson", what = "sp")
medellin <- geojson_read( "otras cuidades/ciudades/medellin.geojson", what = "sp")
cartagena <- geojson_read( "otras cuidades/ciudades/cartagena.geojson", what = "sp")
barranquilla <- geojson_read( "otras cuidades/ciudades/barranquilla.geojson", what = "sp")
bucaramanga <- geojson_read( "otras cuidades/ciudades/bucaramanga.geojson", what = "sp")


# generacion hexahogonos de Cali
cali_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% spTransform(cali, .)
hex_grid <- make_grid(cali_utm, cell_area = 0.1, clip = TRUE)
hex_grid_a<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")%>% spTransform(hex_grid, .)
hex_grid_a$ID<-1:length(hex_grid_a)
geojson_cali<-as.geojson(hex_grid_a)
write(geojson_cali, "otras cuidades/ciudades/cali_subdivi_100.geojson")

# generacion hexahogonos de Medellin
medellin_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% spTransform(medellin, .)
hex_grid <- make_grid(medellin_utm, cell_area = 0.1, clip = TRUE)
hex_grid_a<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")%>% spTransform(hex_grid, .)
hex_grid_a$ID<-1:length(hex_grid_a)
geojson_medellin<-as.geojson(hex_grid_a)
write(geojson_medellin, "otras cuidades/ciudades/medellin_subdivi_100.geojson")


# generacion hexahogonos de Barranquilla
barranquilla_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% spTransform(barranquilla, .)
hex_grid <- make_grid(barranquilla_utm, cell_area = 0.1, clip = TRUE)
hex_grid_a<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")%>% spTransform(hex_grid, .)
hex_grid_a$ID<-1:length(hex_grid_a)
geojson_barranquilla <-as.geojson(hex_grid_a)
write(geojson_barranquilla, "otras cuidades/ciudades/barranquilla_subdivi_100.geojson")

# generacion hexahogonos de Bucaramanga
bucaramanga_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% spTransform(bucaramanga, .)
hex_grid <- make_grid(bucaramanga_utm, cell_area = 0.1, clip = TRUE)
hex_grid_a<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")%>% spTransform(hex_grid, .)
hex_grid_a$ID<-1:length(hex_grid_a)
geojson_bucaramanga<-as.geojson(hex_grid_a)
write(geojson_bucaramanga, "otras cuidades/ciudades/bucaramanga_subdivi_100.geojson")



###########################################################

# función para quitar los hole

setGeneric('dropHole',def=function(poly, ...){
  standardGeneric('dropHole')
})

setMethod('dropHole',signature=signature('Polygon'),
          def=function(poly) {
            #return only Polygons which are not holes
            if (poly@hole) NULL else poly
          }
)

# Drop holey sp::Polygon entries in the @Polygons list of an sp::Polygons class
setMethod('dropHole', signature = signature('Polygons'),
          def = function(poly) {
            noHoles <- lapply(poly@Polygons, dropHole)
            #Remove the NULL entries from the list
            noHoles <- Filter(Negate(is.null), noHoles)
            # Turn back into a (single) Polygons
            # The generator function (sp::Polygons) fills in the other slots!
            # return the new sp:Polygons object
            sp::Polygons(noHoles, ID = poly@ID)
          }
)

# Drop holey parts of sp::Polygons in the @polygons list 
# of an sp::SpatialPolygonsDataFrame
setMethod('dropHole', signature = signature('SpatialPolygonsDataFrame'),
          def = function(poly) {
            noHoles <- lapply(poly@polygons, dropHole)
            # Put the un holey Polygons list back into the @polygons slot 
            poly@polygons <- noHoles
            #return the modified SpatialPolygonsDataFrame 
            poly
          }
)


###########################################################

# generacion hexahogonos de Cartagena

len_hex_grid_a <- 0
cartagena <- geojson_read( "otras cuidades/ciudades/cartagena.geojson", what = "sp")
rm(ab)
# generacion hexahogonos de Cartagena
for(i in 1:length(cartagena@polygons[[1]]@Polygons)){
cat("i:",i,"\n")
cartagena <- geojson_read( "otras cuidades/ciudades/cartagena.geojson", what = "sp")

n <- length(cartagena@polygons[[1]]@Polygons)
vec <- c(1:n);vec <- setdiff(vec,i)

#cartagena@polygons[[1]]@Polygons[c(1:49)] <- NULL 
cartagena@polygons[[1]]@Polygons[vec] <- NULL 
cartagena@polygons[[1]]@plotOrder <- cartagena@polygons[[1]]@plotOrder[3]
cartagena@polygons[[1]]@labpt <- cartagena@polygons[[1]]@Polygons[[1]]@labpt
cartagena@polygons[[1]]@area <- cartagena@polygons[[1]]@Polygons[[1]]@area
cartagena@polygons[[1]]@Polygons[[1]]@hole <- FALSE

cartagena <- dropHole(cartagena)


# generacion hexahogonos de Cartagena

cartagena_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% spTransform(cartagena, .)
hex_grid <- make_grid(cartagena_utm, cell_area = 0.1, clip = TRUE)
hex_grid_a<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")%>% spTransform(hex_grid, .)
hex_grid_a$ID<-(len_hex_grid_a+1):(length(hex_grid_a)+len_hex_grid_a)
len_hex_grid_a <- len_hex_grid_a+length(hex_grid_a)

if(i==1){ab <- hex_grid_a}else{ab<-union(ab,hex_grid_a)}
}

mapview(ab)
geojson_cartagena <-as.geojson(ab)
write(geojson_cartagena, "otras cuidades/ciudades/cartagena_subdivi_100.geojson")

###########################################################

# visualización de los poligonos

cali <- geojson_read( "otras cuidades/ciudades/cali_subdivi_100.geojson", what = "sp")
mapview(cali)
length(cali)


medellin <- geojson_read( "otras cuidades/ciudades/medellin_subdivi_100.geojson", what = "sp")
mapview(medellin)
length(medellin)


barranquilla <- geojson_read( "otras cuidades/ciudades/barranquilla_subdivi_100.geojson", what = "sp")
mapview(barranquilla)
length(barranquilla)


cartagena <- geojson_read( "otras cuidades/ciudades/cartagena_subdivi_100.geojson", what = "sp")
mapview(cartagena)
length(cartagena)
