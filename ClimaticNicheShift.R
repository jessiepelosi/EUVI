######################################
##                                  ##
## Climatic Niche Shifts in EUVI    ##
## Jessie Pelosi, 2025              ##
##                                  ##
######################################    

source("niche.overlap.functions.R")
source("occ.prep.functions.R")
library(ecospat)
library(BIOMOD)
library(ade4)
library(adehabitatHR)
library(sp)
library(gam)
library(MASS)
library(mvtnorm)
library(gbm)
library(dismo)
library(dplyr)
library(sf)

#################################################################################################
############################## preparation of datasets ##########################################
#################################################################################################

# In this section, we prepare the datasets to run the niche analyses

# FIRST: 
# Read-in bioclim layers (each variable is in it's own TIFF, need to read in as a list and then stack)
clim_list <- list.files("./climate/wc2.1_30s/",pattern=".tif",full.names = T)
bio_nums <- as.numeric(gsub(".*_(\\d+)\\.tif$", "\\1", clim_list))  # Extract numeric parts
clim_list <- clim_list[order(bio_nums)]  # Sort by bio number

# Filter the clim list to only include the 5 bioclim variables we are interested in: 
clim_list <- c("./climate/wc2.1_30s//wc2.1_30s_bio_2.tif", "./climate/wc2.1_30s//wc2.1_30s_bio_5.tif", 
               "./climate/wc2.1_30s//wc2.1_30s_bio_6.tif", "./climate/wc2.1_30s//wc2.1_30s_bio_17.tif", 
               "./climate/wc2.1_30s//wc2.1_30s_bio_19.tif")

# stacking the bioclim variables to process them at one go
clim <- raster::stack(clim_list)

# read in the shapefile that contains administrative areas including countries, provinces, and states 
# you will need to use the 'here' library for the following commands 
library(here)
countries <- read_sf("../02_NonNativeFerns/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")

# Subset the shapefile to only include verified areas in the introduced range
EUVI.WUS_shp <- countries %>% 
  filter(name == "California" | name == "Washington" | name == "Oregon" | name == "Nevada" | name == "Idaho")

# Plot the updated shapefile to make sure it looks good 
plot(EUVI.WUS_shp)

# For some reason, using the object that R creates doesn't work for the 'mask' command below
# so we will write out the shapefile and then read back in with the 'shapefile' command 
st_write(EUVI.WUS_shp, 'EUVI.WUS.shp', append = F)
EUVI.WUS_shp <- shapefile("EUVI.WUS.shp")

# This steps crops the stacked climate raster based on the shapefile generated above, 
#but takes a long time to run, so try to only run once. 
clim.WUS <- mask(clim, EUVI.WUS_shp)

# as.data.frame doesn't give the long/lat, need to use rasterToPoints and convert that to df, 
# converting using rasterToPoints uses a lot of storage space, so do each one of these at a time and 
# delete the temp files in between (alternatively run somewhere where you have lots of extra storage space to spare)

clim.WUS_pts <- rasterToPoints(clim.WUS) 
clim.WUS_pts.df <- as.data.frame(clim.WUS_pts) 

# rename the df and columns to be consistent with the remainder of the script
clim1 <- clim.WUS_pts.df
colnames(clim1) <- c("x", "y", "X2", "X5", "X6", "X17", "X19")

# Write out csv so we don't have to do this again! 
write.csv(file = "EUVI.WUS_clim1.csv", clim1, quote = F)

# clim1 is DONE

# SECOND: 

# Subset the shapefile to only include verified areas in the native range
EUVI.EU_shp <- countries %>% 
  filter(admin == "Greece" | admin == "Turkey" | admin == "Bulgaria" | admin == "Iran" | admin == "Georgia" |
           admin == "Ukraine" | admin == "Moldova" | woe_name == "Crimea") 

# Plot the updated shapefile to make sure it looks good 
plot(EUVI.EU_shp)

# For some reason, using the object that R creates doesn't work for the 'mask' command below
# so we will write out the shapefile and then read back in with the 'shapefile' command 
st_write(EUVI.EU_shp, 'EUVI.EU.shp')
EUVI.EU_shp <- shapefile("EUVI.EU.shp")

clim_EU <- mask(clim, EUVI.EU_shp)

# as.data.frame doesn't give the long/lat, need to use rasterToPoints and convert that to df, 
# converting using rasterToPoints uses a lot of storage space, so do each one of these at a time and 
# delete the temp files in between (alternatively run somewhere where you have lots of extra storage space to spare)
clim_EU_pts <- rasterToPoints(clim_EU)
clim_EU_pts.df <- as.data.frame(clim_EU_pts) 

# rename the df and columns to be consistent with the remainder of the script
clim2 <- clim_EU_pts.df
colnames(clim2) <- c("x", "y", "X2", "X5", "X6", "X17", "X19")

# Write out csv so we don't have to do this again! 
write.csv(file = "EUVI.EU_clim2.csv", clim2, quote = F)

# clim2 is DONE 

# THIRD: 

# load climate variable for all site of the study area 1 (column names should be x,y,X1,X2,...,Xn)
clim1<-na.exclude(read.delim("EUVI.WUS_clim1.csv",h=T,sep=","))
clim1 <- dplyr::select(clim1, x, y, X2, X5, X6, X17, X19)

# load climate variable for all site of the study area 2 (column names should be x,y,X1,X2,...,Xn)
clim2<-na.exclude(read.delim("EUVI.EU_clim2.csv",h=T,sep=",")) 
clim2 <- dplyr::select(clim2, x, y, X2, X5, X6, X17, X19)


# global climate for both ranges
clim12<-rbind(clim1,clim2)
############################################
## Occurrence data filtering and cleaning ## 
############################################

# Note that most cleaning was performed by hand by B.S. Barker, for the native range
EUVI.EU.occ <- read.delim("EUVI_Eurasia_records.tsv") 
EUVI.EU.occ <- EUVI.EU.occ %>% 
  filter(Notes != "Record for province, not precise")

# Add in any occ data from GBIF: 
library(rgbif)
EUVI.gbif <- occ_data(taxonKey = "4465015")
EUVI.gbif.data <- EUVI.gbif$data 

## Get occurrence data from the sampling in this study: 
env <- read.delim("individual_climate_data_30s.txt", sep = " ")
env$Site <- as.character(env$Site)

# Filter for individuals/pops in the GEA dataset 
locations <- read.delim("EUVI_locations.txt", header =T)
# Filter locations to only samples that we are using: 
popmap3 <- read.delim("samples_to_keep.txt", sep = " ")

locations_update <- distinct(locations[locations$EUVI_Sample %in% popmap3$V1,])

# Add locations from B.S.B. on WUS releases: 
WUS_releases <- read.csv("EUVI_WUS_releases - Sheet1.csv")
WUS_releases <- WUS_releases %>% 
  filter(County_level == 0)
WUS_releases_coords <- dplyr::select(WUS_releases, Latitude, Longitude)
colnames(WUS_releases_coords) <- c("Lat", "Long")

# Filter to all samples we are keeping: 
env <- env[env$EUVI_Sample %in% locations_update$EUVI_Sample,]
WUS_releases_coords$Lat <- as.numeric(WUS_releases_coords$Lat)
WUS_occs <- full_join(env, WUS_releases_coords)


# Now filter to only include samples in the introduced range and only the variables we are interested in
env_WUS <- WUS_occs %>% 
  dplyr::filter(Long < 0) %>% 
  dplyr::select(Site, Lat, Long, bio5, bio6, bio19, bio17, bio2) %>% 
  distinct()

# Now let's merge all these dfs together. 
# First, need to trim to only the lat/long, and rename so all the columns are the same across dfs so we can merge
EUVI.EU.occ.trim <- EUVI.EU.occ %>% 
  mutate(x=Longitude, y= Latitude) %>% 
  dplyr::select(x,y) %>% 
  na.omit()

EUVI.gbif.data.trim <- EUVI.gbif.data %>% 
  mutate(x=decimalLongitude, y=decimalLatitude) %>% 
  dplyr::select(x,y) %>% 
  na.omit()

EUVI.env.WUS.trim <- env_WUS %>% 
  mutate(x=Long, y = Lat) %>% 
  dplyr::select(x,y) %>% 
  na.omit()

EUVI.EU.occ.trim$y <- as.double(EUVI.EU.occ.trim$y)
EUVI.gbif.data.trim$y <- as.double(EUVI.gbif.data.trim$y)

tmp <- rbind(EUVI.EU.occ.trim, EUVI.gbif.data.trim)

all.EUVI.occs <- rbind(tmp, EUVI.env.WUS.trim)
write.csv(file = "EUVI_occs.csv", quote = F, row.names = F, all.EUVI.occs)


# loading occurrence sites for the species (column names should be x,y)
# we want to use 'select' first and then 'na.exclude' otherwise we drop all the rows b/c there are other rows that have NAs
occ <-read.delim("EUVI_occs.csv",h=T,sep=",")

# These lines will estimate the Mantel autocorrelation and the optimal distance to avoid spatial autocorrelation
# which we will use to thin the occurrence data
ecospat.mantel.correlogram(dfvar = clim1, colxy = c("x","y"), n = 100, colvar = c("X2","X5","X6","X17","X19"), max = 1000, nclass=10, nperm=100)
ecospat.mantel.correlogram(dfvar = clim2, colxy = c("x","y"), n = 100, colvar = c("X2","X5","X6","X17","X19"), max = 1000, nclass=10, nperm=100)

# From these graphs, it looks like 200m is the minimum distance needed to thin the occurrence data to avoid spatial autocorrelation 


# remove occurrences closer than a minimum distance to each other (remove aggregation). 
# Setting min.dist=0 will remove no occurrences.
# From what I can tell, this function is based on Euclidean distances of lat and long in decimals,
# So we need to convert our minimum distance from above in meters first to something equivalent in lats and longs 

# We know that one degree of latitude is equal to 111.111 km, so to go from meters to degrees of latitutde
# We should divide distance above (0.2km) by 111.111km by the the 
# We need to know the lat to calculate the long, so let's simply use the latitude value as our minimum distance in degrees for the next function

lat <- 0.2/111.111

occ.sp<-occ.desaggragation(df=occ,colxy=1:2,min.dist=lat,plot=F)

# Try using spThin -- this does the same thing as the line above, I just wanted to check that I was doing the distance right
#library(spThin)
#occ$species <- "EUVI"
#thin(occ, lat.col = "x", long.col = "y", spec.col = "species", thin.par = 0.2, reps = 100)

# create sp occurrence dataset by adding climate variables from the global climate datasets
# resolution should be the resolution of the climate data grid

# Check that resolution is the right one (matches resolution of the climate data grid) 
occ.sp1<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim1,colvarxy=1:2,colvar="all",resolution=0.16666))
occ.sp2<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim2,colvarxy=1:2,colvar="all",resolution=0.16666))


#################################################################################################
################### row weighting and grouping factors for ade4 functions  ######################
#################################################################################################

row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))
# Here, I write out the dataframe to csv in case of any problems, so I can re-load if necessary rather than re-run. 
data.table::fwrite(as.data.frame(row.w.occ), "row.w.occ.csv")

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
# Here, I write out the dataframe to csv in case of any problems, so I can re-load if necessary rather than re-run. 
data.table::fwrite(as.data.frame(row.w.env), "row.w.env.csv")

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

# Only need to do these next steps once
# rbind takes a lot of RAM, so this is where I switched to using the HPC
# I also use 'fwrite' here from the 'data.table' package to write out the csv 
# since it is more memory-efficient than 'write.csv'
Xvar<-c(3:7)
data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
head(data.env.occ)
# Here, I write out the dataframe to csv in case of any problems, so I can re-load if necessary rather than re-run. 
data.table::fwrite(data.env.occ, "data.env.occ.csv")

row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))


pca.env <- ade4::dudi.pca(data.env.occ, row.w = row.w.env, scannf=F, nf=2, center = T, scale = T)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
pca.env$co

# PCA scores for the species native distribution
scores.globclim<- pca.env$li[row.clim12,]
# PCA scores for the species invasive distribution
scores.clim.inv<- pca.env$li[row.clim1,]
# PCA scores for the species native distribution
scores.clim.nat<- pca.env$li[row.clim2,]
#PCA scores for the species occurrence points
scores.sp.inv<- pca.env$li[row.sp1,]
scores.sp.nat<- pca.env$li[row.sp2,]

grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.nat, sp=scores.sp.nat, R=100, th.sp = 0)
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.inv, sp=scores.sp.inv, R=100, th.sp = 0)

# D reports the % of overlap between the native and invaded niches
D.overlap <- ecospat.niche.overlap(grid.clim.nat, grid.clim.inv, cor =T)$D
D.overlap

#"Niche equivalency test H1: the observed overlap between the native and invaded niche is higher than if the two niches are randomized, 
#the niche expansion of the invaded niche is lower, the niche stability is higher and the niche unfilling is lower than if the two niches are randomized."

#eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep=1000, intersection = 0.1, 
#                                          overlap.alternative = "higher", 
#                                          expansion.alternative = "lower", 
#                                          stability.alternative = "higher", 
#                                          unfilling.alternative = "lower")

#ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

# Perform niche similarity test, rec'd to do 1000 reps

# "Niche similarity test H1: the observed overlap between the native and invaded is higher than randomly shifted invasive niches in the 
# invaded study area, the niche expansion of the invaded niche is lower, the niche stability is higher and the niche unfilling is lower than if the two niches are randomized."

sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=1000, overlap.alternative = "lower", 
                                          expansion.alternative = "higher",
                                          stability.alternative = "lower",
                                          unfilling.alternative = "higher",
                                          intersection = 0.1,
                                          rand.type=1)

ecospat.plot.overlap.test(sim.test, "D", "Similarity")


ecospat.plot.overlap.test(sim.test, "D", "Similarity")
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")

ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv, intersection =0.25)

niche.dyn <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv)

ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant = 0.25, interest = 2,
                       name.axis1 = "PC1", name.axis2 = "PC2", 
                       col.unf = "deepskyblue3", col.stab = "firebrick4", col.exp = "orange", col.abn = "skyblue1",
                       colZ1 = "deepskyblue3", colZ2 = "orange", transparency = 50)

library(dplyr)
library(ggplot2)
source("./functions.R")

grid.clim.pooled <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.globclim, 
                                          sp=rbind(scores.sp.nat, scores.sp.inv),
                                          R=100, th.sp=0)

# This will give different values b/c we are dealing with analog and non-analog space
calculate.aban.un.stab.exp.pion(grid.clim.pooled = grid.clim.pooled, grid.clim.nat = grid.clim.nat, 
                                grid.clim.inv = grid.clim.inv, remove.marginal = 0.25)

# let's try to re-plot with ggplot using Riera's code:
cat <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv, intersection = 0.25)$dyn
cat.df <- terra::as.data.frame(cat, xy=T) %>% mutate(fill.ausep = case_when(lyr.1 == 0 ~"white", 
                                                                            lyr.1 == 1 ~ "skyblue1",
                                                                            lyr.1 == 2 ~ "deepskyblue3",
                                                                            lyr.1 == 3 ~ "firebrick4",
                                                                            lyr.1 == 4 ~ "orange",
                                                                            lyr.1 == 5 ~ "gold2",
                                                                            lyr.1 == 6 ~ "gray"))

grid.clim.inv$z.uncor[grid.clim.inv$z.uncor[] == 0] <- NA
grid.inv.z.uncor.df <- terra::as.data.frame(grid.clim.inv$z.uncor, xy = T)

grid.invaded.Z.df <- terra::as.data.frame(grid.clim.inv$Z, xy = T)
grid.native.Z.df <- terra::as.data.frame(grid.clim.nat$Z, xy = T)

plot.ssp <- ggplot() +
  geom_raster(data = cat.df, aes(x = x, y = y, fill = fill.ausep), show.legend = F) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  geom_raster(data = grid.inv.z.uncor.df, aes(x = x, y = y, fill = lyr.1), alpha = 0.5, show.legend = F) +
  scale_fill_gradient(low= "white", high = "black") +
  # add contours of the niche in invaded range
  geom_contour(data = grid.invaded.Z.df, aes(x = x, y = y, z = lyr.1), breaks = quantile(grid.clim.inv$Z[grid.clim.inv$Z > 0], probs = 0), color = "orange", linetype = "solid") +
  geom_contour(data = grid.invaded.Z.df, aes(x = x, y = y, z = lyr.1), breaks = quantile(grid.clim.inv$Z[grid.clim.inv$Z > 0], probs = 0.25), color = "orange", linetype = "dashed") +
  # add contours of the niche in the native range
  geom_contour(data = grid.native.Z.df, aes(x = x, y = y, z = lyr.1), breaks = quantile(grid.clim.nat$Z[grid.clim.nat$Z > 0],probs = 0), color = "deepskyblue3", linetype = "solid") +
  geom_contour(data = grid.native.Z.df, aes(x = x, y = y, z = lyr.1), breaks = quantile(grid.clim.nat$Z[grid.clim.nat$Z > 0],probs=0.25), color = "deepskyblue3", linetype = "dashed") +
  #geom_segment(aes(x = median(scores.sp.nat$Axis1), y = median(scores.sp.nat$Axis2), xend = median(scores.sp.inv$Axis1), yend = median(scores.sp.inv$Axis2), arrow = arrow(length=unit(0.1, "cm")), linewidth = 1)) +
  #geom_segment(aes(x = median(scores.clim.nat$Axis1), y = median(scores.clim.nat$Axis2), xend = median(scores.clim.inv$Axis1), yend = median(scores.clim.inv$Axis2), arrow = arrow(length=unit(0.1, "cm")), linewidth = 1, color = "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand =c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "PC1 (51.64% of variation)", y = "PC2 (26.44 % of variation)", title = "") #+
#  theme(plot.title = element_text(face = "italic", size = 14))

plot.ssp

ggsave("EUVI_ECOSPAT.pdf", height = 5, width = 7)


#niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv)
#ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest= 2, title= "Niche Overlap", name.axis1="PC1",
#                       name.axis2="PC2", 
#                       col.unf = "skyblue2", col.stab = "firebrick4", col.exp = "orange", 
#                       colZ1 = "skyblue2", colZ2 = "orange", transparency =52)

#ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv, col = "black")

#ecospat.plot.overlap.test(sim.test, "D", "Similarity")
#ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
#ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
#ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")

