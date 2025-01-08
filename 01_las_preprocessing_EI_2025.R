rm(list=ls(all=TRUE))
library(lidR)
library(future)
library(RCSF)
library(terra)
library(lidRmetrics)
library(sf)
library(ForestTools)
library(raster)
library(rgdal)

setwd("G:/Project2023/02_LidR")
outpath = "G:/Project2023/02_LidR" 
LASfile="G:/Project2023/01_LiDAR/A_01.las"

plan(multisession, workers=7L)
ctg = readLAScatalog(LASfile)
plot(ctg)

opt_chunk_buffer(ctg) <- 10
opt_chunk_size(ctg) <- 250

opt_output_files(ctg) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_tiled")
newctg <- catalog_retile(ctg)
plot(newctg,chunk=TRUE)

# ---- Clip ----
shp <- st_read("G:/Project2023/00_Shapefile/ROI.shp")
opt_output_files(newctg) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_roi")
cropped_las <- clip_roi(newctg, shp)

# ---- density calculation 1----
grid_density(cropped_las, res = 1)
density1 <- list.files(path = outpath, pattern = '_roi.tif$', full.names = T)
density1_mosaic <- vrt(density1, overwrite = T)
plot(density1_mosaic)

mean <- mean(density1_mosaic[], na.rm = TRUE)
sd <- sd(density1_mosaic[], na.rm = TRUE)
maximum<- mean + sd

# ---- thin ----
opt_output_files(cropped_las) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_thinned")
thinned <- decimate_points(cropped_las, homogenize(mean,1))

# ---- density calculation 2----
grid_density(thinned, res = 1)
density2 <- list.files(path = outpath, pattern = '_thinned.tif$', full.names = T)
density2_mosaic <- vrt(density2, overwrite = T)
plot(density2_mosaic)

# ---- Classify noise points ----
opt_output_files(thinned) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_classifynoise")
classified_ctg1<- classify_noise(thinned, ivf(5,6)) 

# ---- Classify ground points ----
opt_output_files(classified_ctg1) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_classifyground")
classified_ctg2 <- classify_ground(classified_ctg1, csf(class_threshold = 0.5, cloth_resolution = 0.1, rigidness = 2)) 
plot(classified_ctg2,chunk=TRUE) #The points classified as 'ground' are assigned a value of 2

# ----- DTM -----
opt_output_files(classified_ctg2) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_dtm")
rasterize_terrain(classified_ctg2, res = 0.1, tin())
dtm_tiles <- list.files(path = outpath, pattern = '_dtm.tif$', full.names = T)
dtm_mosaic <- vrt(dtm_tiles, overwrite = T)
dtm_mosaic
writeRaster(dtm_mosaic, filename = 'dtm_mosaic.tif', overwrite = T)

# ----- DTM smooth-----
dtm_smooth <- dtm_mosaic %>%focal(w = matrix(1, 25, 25), fun = base::mean, na.rm = TRUE,pad = TRUE)
writeRaster(dtm_smooth, filename = 'dtm_mosaic_smooth.tif', overwrite = T)
plot(dtm_mosaic, bg = "white") 
dtm_prod <- terra::terrain(dtm_mosaic, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col = gray(0:50/50), legend = FALSE)
plot(dtm_smooth, bg = "white") 
dtm_prod_smooth <- terra::terrain(dtm_smooth, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade_smooth <- terra::shade(slope = dtm_prod_smooth$slope, aspect = dtm_prod_smooth$aspect)
plot(dtm_hillshade_smooth, col = gray(0:50/50), legend = FALSE)

# ---- Normalize point cloud ----
opt_output_files(classified_ctg2) <-  paste0(outpath, "/{XLEFT}_{YBOTTOM}_norm_tin")
opt_filter(classified_ctg2) <- '-drop_withheld'
opt_filter(ctg) <- "-drop_class 18"
ctg_norm_tin <- normalize_height(classified_ctg2, tin())

# ----- Canopy Height Model -----
opt_output_files(ctg_norm_tin) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_chm")
rasterize_canopy(ctg_norm_tin, res = 0.1,algorithm = p2r(na.fill = knnidw()))
chm_tiles <- list.files(path = outpath, pattern = '_chm.tif$', full.names = T)
chm_mosaic <- vrt(chm_tiles, overwrite = T)

#--- Remove extreme CHM values ---
chm_mosaic <- clamp(chm_mosaic, 0, 50, values = TRUE)
#--- Set layer name to Z ---
names(chm_mosaic) <- 'Z'
# ----- Fill CHM -----
chm_filled <- focal(chm_mosaic,w = 3,fun = "mean",na.policy = "only",na.rm = TRUE)
names(chm_filled) <- 'Z'
# ----- Smooth CHM -----
fgauss <- function(sigma, n) {
  m <- matrix(ncol = n, nrow = n)
  col <- rep(1:n, n)
  row <- rep(1:n, each = n)
  x <- col - ceiling(n/2)
  y <- row - ceiling(n/2)
  m[cbind(row, col)] <- 1/(2 * pi * sigma^2) * exp(-(x^2 + y^2)/(2 * sigma^2))
  m/sum(m)
}
chm_smooth <- terra::focal(chm_filled,w = fgauss(1, n = 5))
names(chm_smooth) <- 'Z'
plot(chm_mosaic, col = height.colors(50))
plot(chm_filled, col = height.colors(50))
plot(chm_smooth, col = height.colors(50))
writeRaster(chm_mosaic, filename = 'chm_mosaic.tif', overwrite = T)
writeRaster(chm_filled, filename = 'chm_filled.tif', overwrite = T)
writeRaster(chm_smooth, filename = 'chm_smooth.tif', overwrite = T)

# ----- Tree Crown Delineation ----- 
image.fn <- "chm_smooth.tif"
CHM_input <- raster(image.fn)
shapefile <- "G:/Project2023/00_Shapefile/ROI.shp"
shp <- st_read(shapefile)
CHM <- crop(CHM_input, shp)
plot(CHM, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
lin <- function(x){x * 0.05 + 0.6}
ttops <- vwf(CHM, winFun = lin, minHeight = 1.5)
plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)
crownsPoly <- mcws(treetops = ttops, CHM, format = "polygons", minHeight = 1.5, verbose = FALSE)
plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)
writeOGR(crownsPoly, "G:\\Project2023\\00_Shapefile", "Crown", driver = "ESRI Shapefile")

# ----- Extract Invidual Tree 
shapefile <- st_read("G:/Project2023/00_Shapefile/Crown_modified.shp")
shapefile_no_z <- st_zm(shapefile)
opt_output_files(ctg_norm_tin) <- paste0(outpath, "/{XLEFT}_{YBOTTOM}_Tree")
opt_filter(ctg_norm_tin) <- "-drop_withheld -drop_z_below 0"
opt_filter(ctg_norm_tin) <- "-drop_class 18 2"
cropped_las <- clip_roi(ctg_norm_tin, shapefile_no_z)

# ----- Calculate Lidar Metrics
las_dir <- "G:\\Project2023\\02_LidR"
las_files <- list.files(las_dir, full.names = TRUE, pattern = "\\_Tree.las$", recursive = FALSE)
results_list <- list()
cat("Processing", length(las_files), "files...\n")

# Initialize counter for unique tree IDs
tree_counter <- 1

for (las_file_path in las_files) {
  cat("Processing:", basename(las_file_path), "\n")
  las <- readLAS(las_file_path)
  las@data$treeID <- rep(tree_counter, nrow(las))  
  metrics <- tree_metrics(las, ~lidRmetrics::metrics_set3(X, Y, Z, i = Intensity, 
                                                          ReturnNumber = ReturnNumber, 
                                                          NumberOfReturns = NumberOfReturns,
                                                          vox_size = 0.1))
  results_list[[las_file_path]] <- metrics
  tree_counter <- tree_counter + 1  
}

combined_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  df <- as.data.frame(results_list[[i]])
  df$Element_ID <- i  
  return(df)
}))

write.csv(combined_df, "combined_data_0.1.csv", row.names = FALSE)

