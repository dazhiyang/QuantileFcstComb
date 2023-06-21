#################################################################################
# This code is written by Dazhi Yang
# Department of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "SolarData", "ncdf4", "camsRad", "doSNOW", "zoo")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir.ecmwf.c <- "/Volumes/Data/ECMWF_Data_Article/ECMWF/12Z" # dir for ecmwf control forecast
dir.ecmwf.p <- "/Volumes/Data/ECMWF ENS 12Z/12Z" # dir for ecmwf perturbed forecast
dir.surfrad <- "/Volumes/Macintosh Research/Data/SURFRAD" # dir for SURFRAD data
dir.nsrdb <- "/Users/dyang/Dropbox/Working papers/Ensemble/Data/NSRDB" # dir for NSRDB data
dir.mcclear <- "/Users/dyang/Dropbox/Working papers/Ensemble/Data/McClear" # dir for McClear data
dir.processed <- "/Users/dyang/Dropbox/Working papers/Ensemble/Data/Processed" # dir for processed data
cams_set_user("yangdazhi.nus@gmail.com")
yr <- c(2017:2020)
stn <- "tbl" # process on station at a time, bon, dra, fpk, gwn, psu, sxf, tbl
#################################################################################

# read metadata of SURFRAD stations
loc <- SURFRAD.loc
# index of the station in the "loc" data.frame
stn.ind <- match(stn, loc$stn)
# time zone of the station
tz <- loc$tz[stn.ind]

#################################################################################
# get PSM3 data from NSRDB (only run once)
################################################################################
# # for each year, download the hourly files
# for(i in 1:length(yr))
# {
#   PSM.get(lon = loc$lon[stn.ind],
#           lat = loc$lat[stn.ind],
#           api.key = "8JRR5hfacHecdknYn5lE7BoCl3cci80IfaPGB2WC",
#           attributes = "ghi,clearsky_ghi,solar_zenith_angle",
#           name = "Dazhi+Yang",
#           affiliation = "HIT",
#           year = as.character(yr[i]),
#           leap.year = "true",
#           interval = "60",
#           utc = "true",
#           reason.for.use = "research",
#           email = "yangdazhi.nus@gmail.com",
#           mailing.list = "false",
#           directory = file.path(dir.nsrdb, stn))
#   cat("download for year", yr[i], "at", stn,  "done \n") # print progress on the screen
#   Sys.sleep(1) # sleep 1 sec, as there is a rate limit for downloading NSRDB
# }

setwd(file.path(dir.nsrdb, stn)) # set the directory to NSRDB, station
files <- dir(pattern = "*.csv") # get files
files <- files[which(substr(files, nchar(files)-7, nchar(files)-4) %in% yr)]

# use "do.call" to read all csv files, and concatenate them
NSRDB <- do.call("rbind", lapply(files, read.csv, header = TRUE, skip = 2, sep = ","))

# convert the data.frame into tibble
# join different time elements and make DataTime object
# select only the relevant variables
# change the variable names
# shift the time for 30 min, such that the time stamp aligns with that of SURFRAD
NSRDB <- as_tibble(NSRDB) %>%
  dplyr::mutate(Time = lubridate::ymd_hm(paste(paste(Year, Month, Day, sep = "-"), paste(Hour, Minute, sep = ":"), sep = " "))) %>%
  dplyr::select(c("Time", "GHI", "Clearsky.GHI", "Solar.Zenith.Angle")) %>%
  dplyr::rename(Ys = GHI, REST2 = Clearsky.GHI, Z = Solar.Zenith.Angle) %>%
  mutate(Time = Time + 30*60)
# four columns: Time, satellite-derived GHI (sghi), clear-sky satellite-derived GHI (sghic). and zenith angle (szen)
NSRDB 

#################################################################################
# get McClear clear-sky irradiance
#################################################################################
# # use camsRad package to get McClear
# cams_set_user("yangdazhi.nus@gmail.com")
# McClear <- cams_get_mcclear(lat = loc$lat[stn.ind],
#                             lng = loc$lon[stn.ind],
#                             date_begin = paste(yr[1],"-01-01", sep = ""),
#                             date_end = paste(yr[length(yr)],"-12-31", sep = ""),
#                             time_step = "PT15M",
#                             alt = loc$elev[stn.ind],
#                             verbose = FALSE)
# 
# # put McClear data frame into a tibble
# McClear <- as_tibble(McClear)
# # arrange the dataset, takes only the clear-sky GHI values
# # because of 15-min data, times 4 to GHI, to get W/m2
# McClear <- McClear %>% mutate(Time = ymd_hms(timestamp)) %>%
#   mutate(McClear = `Clear sky GHI`) %>% dplyr::select(-c(1:2, 3:6)) %>%
#   mutate(McClear = McClear*4) %>%
#   mutate(McClear = round(McClear, 3))
# # save to disk
# setwd(file.path(dir.mcclear, stn)) # set the directory to McClear, station
# write.table(McClear, file = paste0(loc$stn[stn.ind], "_15min.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# read McClear and aggregate to 1 h
setwd(file.path(dir.mcclear, stn)) # set the directory to McClear, station
file <- dir(pattern = "*.txt") # get files
McClear <- tibble(read.table(file, sep = "\t", header = TRUE))
McClear <- McClear %>%
  mutate(Time = ceiling_date(ymd_hms(Time), "1 h")) %>%
  group_by(Time) %>%
  summarise_all(., mean) %>%
  ungroup()

#################################################################################
# read and aggregate SURFRAD data
#################################################################################
# initialize empty object to hold data from one station
surfrad <- NULL
for(i in 1:length(yr)) 
{
  # SURFRAD data is in daily txt files
  # set directory and get the file names
  setwd(file.path(dir.surfrad, stn, yr[i]))
  files.tmp <- dir(pattern = ".dat")
  
  # read and pre-process SURFRAD data
  tmp <- SURFRAD.read(files = files.tmp, 
                      directory = file.path(dir.surfrad, stn, yr[i]), 
                      use.original.qc = FALSE, 
                      use.qc = TRUE, 
                      test = c("ext", "closr", "dr"), 
                      progress.bar = TRUE, 
                      agg = 60)
  # due to aggregation, there would be one repeated row, remove that
  # and then, combine with the main data.frame
  surfrad <- rbind(surfrad, tmp[-nrow(tmp),])
  
  print(i)
}

# select relevant variables
# rename dw_solar as ghi
# round the data to 3 decimal digits
surfrad <- surfrad %>%
  dplyr::select(one_of("Time", "dw_solar", "Ioh")) %>%
  dplyr::rename(Yg = dw_solar, E0 = Ioh) %>%
  mutate_if(., is.numeric, round, 3)
# four columns: Time, zenith angle (zen), ghi (ghi), extraterrestrial irradiance (Ioh)
surfrad

#################################################################################
# Get ECMWF-HRES (i.e., control forecasts)
#################################################################################
# find out the collocated pixel
setwd(dir.ecmwf.c)
files <- dir(pattern = "*.nc", recursive = TRUE)
ncin <- nc_open(files[1]) # open nc
collocate.lon.index <- which.min(abs(loc$lon[stn.ind]-ncvar_get(ncin, "longitude")))
collocate.lat.index <- which.min(abs(loc$lat[stn.ind]-ncvar_get(ncin, "latitude")))
nc_close(ncin) # close nc

# set up parallel processing
# number of cores used for parallel
cl <- makeCluster(8) 
registerDoSNOW(cl)

#foreach execution for variables
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
CTRL <- foreach(j = 1:length(files), .combine = 'rbind', .options.snow = opts, .packages = c("ncdf4", "dplyr")) %dopar% {
  # open nc
  ncin <- nc_open(files[j]) 
  # get time stamps
  time <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + 3600*ncvar_get(ncin, "time")
  # get variables
  ssrd <- ncvar_get(ncin, "ssrd")[collocate.lon.index, collocate.lat.index,]
  t2m <- ncvar_get(ncin, "t2m")[collocate.lon.index, collocate.lat.index,]
  u10 <- ncvar_get(ncin, "u10")[collocate.lon.index, collocate.lat.index,]
  v10 <- ncvar_get(ncin, "v10")[collocate.lon.index, collocate.lat.index,]
  fal <- ncvar_get(ncin, "fal")[collocate.lon.index, collocate.lat.index,]
  # close nc
  nc_close(ncin) 
  
  #########################
  # unit conversion
  #########################
  # ssrd is in [J/m2], which is also an accumulated parameter, use difference
  ssrd[2:length(ssrd)] <- diff(ssrd)/3600 # [W/m2]
  t2m <- t2m - 273.15 # degC
  w10 <- sqrt(u10^2 + v10^2)
  #########################
  
  # get 12-35-h-ahead forecasts, since index 1 is analysis, start from 13
  time <- time[13:36]
  ssrd <- round(ssrd[13:36],3)
  t2m <- round(t2m[13:36],3)
  w10 <- round(w10[13:36],3) 
  fal <- round(fal[13:36],3) 
  ctrl <- bind_cols(tibble(Time = time), ssrd_C = ssrd, t2m_C = t2m, w10_C = w10, fal_C = fal)
  ctrl
}
close(pb)
stopCluster(cl)


#################################################################################
# Get ECMWF-ENS
#################################################################################
# find out the collocated pixel
setwd(dir.ecmwf.p)
files <- dir(pattern = "*.nc", recursive = TRUE)
ncin <- nc_open(files[1]) # open nc
collocate.lon.index <- which.min(abs(loc$lon[stn.ind]-ncvar_get(ncin, "longitude")))
collocate.lat.index <- which.min(abs(loc$lat[stn.ind]-ncvar_get(ncin, "latitude")))
nc_close(ncin) # close nc

# set up parallel processing
# number of cores used for parallel
cl <- makeCluster(8) 
registerDoSNOW(cl)

#foreach execution for variables
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ENS <- foreach(j = 1:length(files), .combine = 'rbind', .options.snow = opts, .packages = c("ncdf4", "dplyr")) %dopar% {
  # open nc
  ncin <- nc_open(files[j]) 
  # get time stamps
  time <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + 3600*ncvar_get(ncin, "time")
  # get variables
  ssrd <- ncvar_get(ncin, "ssrd")[collocate.lon.index, collocate.lat.index,,]
  ssrd <- t(ssrd)
  colnames(ssrd) <- paste0("ssrd_M", 1:50)
  # close nc
  nc_close(ncin) 
  
  #########################
  # unit conversion
  #########################
  # ssrd is in [J/m2], which is also an accumulated parameter, use difference
  ssrd[2:nrow(ssrd),] <- apply(ssrd, 2, diff)/3600 # [W/m2]
  #########################
  
  # get 12-35-h-ahead forecasts, since index 1 is analysis, start from 13
  time <- time[13:36]
  ssrd <- round(ssrd[13:36, ],3)
  ens <- bind_cols(tibble(Time = time), as_tibble(ssrd))
  ens
}
close(pb)
stopCluster(cl)

#################################################################################
# Combine everything together
#################################################################################
# combine SURFRAD, McClear, and NSRDB using left_join
data <- surfrad %>%
  left_join(., McClear, by = "Time") %>%
  left_join(., NSRDB, by = "Time") %>%
  dplyr::select(one_of("Time", "Yg", "Ys", "REST2", "McClear", "E0", "Z")) 

# combine data with CTRL and ENS forecast
data <- data %>%
  left_join(., CTRL, by = "Time") %>%
  left_join(., ENS, by = "Time")

# write data to file
setwd(dir.processed)
file.name <- paste0(loc$stn[stn.ind],"_ECMWF_", yr[1], "_", yr[length(yr)], ".txt")
write.table(data, file = file.name, quote = FALSE, row.names = FALSE, sep = "\t")

# check alignment
period <- 11000:11100
plot(data$REST2[period], type = "l")
lines(data$McClear[period], col = 2)
lines(data$ssrd_C[period], col = 3)
lines(data$Yg[period], col = 4)
lines(data$Ys[period], col = 5)

plot(data$Yg, data$ssrd_C, pch = ".")
abline(0,1,col=2)
plot(data$Ys, data$ssrd_M1, pch = ".")
abline(0,1,col=2)


