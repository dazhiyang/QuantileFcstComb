#################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
# load the required packages
libs <- c("dplyr", "lubridate", "SolarData", "ensembleMOS", "quantreg", "qrnn", "quantregForest")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/dyang/Dropbox/Working papers/Combine/Data" # dir for processed data
stns <- c("bon", "dra", "fpk", "gwn", "psu", "sxf", "tbl")
yr.tr <- c(2017,2018) # training years
yr.te <- c(2019,2020) # verification years
Q <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99) # a set of quantiles to evaluate
m <- length(Q) # number of ensemble members
# function to calculate pinball loss
pinball <- function(y, f, tau)
{
  pb <- matrix(NA, nrow = nrow(f), ncol = ncol(f))
  for(j in 1:length(tau))
  {
    q <- f[,j]
    pb[,j] <- ifelse(y >= q, tau[j] * (y-q), (1-tau[j]) * (q-y))
  }
  mean(rowMeans(pb))
}
#################################################################################

# read metadata of SURFRAD stations
loc <- SURFRAD.loc
# reorder rows of "loc", to follow the order of the stations
loc <- loc[match(stns, loc$stn),]

#################################################################################
# Make post-processed forecasts
#################################################################################
setwd(dir0)

pb = txtProgressBar(min = 0, max = length(stns), initial = 0, style = 3) 
for(i in 1:length(stns))
{
  # read the arranged dataset as provided in the supplementary materials
  # convert time to lubridate time format
  # change from UTC to local time
  # remove the first and last day, which has incomplete data due to time zone change
  # fill missing SURFRAD observations with NSRDB irradiance
  # remove unphysical negative irradiance forecasts
  # remove nighttime data (those zenith angles > 85)
  # select only relevant columns
  assign(stns[i], value = tibble(read.table(file.path(dir0, paste0(stns[i],"_ECMWF_2017_2020.txt")), sep = "\t", header = TRUE)))
  data <- get(stns[i]) %>%
    mutate(Time = ymd_hms(Time)) %>%
    mutate(Time = Time + loc$tz[i]*3600) %>%
    filter(date(Time) >= date("2017-01-02") & date(Time) <= date("2020-12-30")) %>%
    mutate(Yg = ifelse(is.na(Yg), Ys, Yg)) %>%
    mutate(ssrd_C = ifelse(ssrd_C < 0, 0, ssrd_C)) %>%
    filter(Z < 85) %>%
    dplyr::select(one_of("Time", "Yg", "REST2", names(.)[which(grepl("_M", names(.)))]))
  
  # train/test split 
  data.tr <- data %>% filter(year(Time) %in% yr.tr)
  data.te <- data %>% filter(year(Time) %in% yr.te)
  
  # make empty tibble to hold post-processed forecasts
  # set up ensembleData, which is required by the two packages, namely, "ensembleBMA" and "ensembleMOS"
  train.eD <- ensembleData(forecasts = t(apply(data.tr[,which(colnames(data.tr) %in% paste0("ssrd_M", c(1:50)))], 1, sort)), dates = data.tr$Time, observations = data.tr$Yg, forecastHour = 12, initializationTime = "12", exchangeable = rep(1, 50))
  test.eD <- ensembleData(forecasts = t(apply(data.te[,which(colnames(data.te) %in% paste0("ssrd_M", c(1:50)))], 1, sort)), dates = data.te$Time, observations = data.te$Yg, forecastHour = 12, initializationTime = "12", exchangeable = rep(1, 50))
  
  #####################################
  # post-process forecasts
  #####################################
  # make empty tibble to hold post-processed forecasts
  fcst <- data.te %>% 
    dplyr::select(one_of("Time", "Yg", "REST2")) %>% 
    rename(y = Yg) %>% 
    mutate(y = round(y, 1))
  
  # EMOS with Gaussian CRPS minimization
  EMOS1.fit <- fitMOS(train.eD, model = "normal", control = controlMOSnormal(scoringRule = c("crps")))
  EMOS1 <- quantileForecast(EMOS1.fit, test.eD, quantiles = Q)
  EMOS1 <- round(apply(EMOS1, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  colnames(EMOS1) <- paste0("EMOS1.", Q*100)
  fcst <- bind_cols(fcst, EMOS1)
  
  # EMOS with Gaussian IGN minimization (i.e., ML estimation)
  EMOS2.fit <- fitMOS(train.eD, model = "normal", control = controlMOSnormal(scoringRule = c("log")))
  EMOS2 <- quantileForecast(EMOS2.fit, test.eD, quantiles = Q)
  EMOS2 <- round(apply(EMOS2, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  colnames(EMOS2) <- paste0("EMOS2.", Q*100)
  fcst <- bind_cols(fcst, EMOS2)
  
  # EMOS with lognormal CRPS minimization
  EMOS3.fit <- fitMOS(train.eD, model = "lognormal", control = controlMOSlognormal(scoringRule = c("crps")))
  EMOS3 <- quantileForecast(EMOS3.fit, test.eD, quantiles = Q)
  EMOS3 <- round(apply(EMOS3, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  colnames(EMOS3) <- paste0("EMOS3.", Q*100)
  fcst <- bind_cols(fcst, EMOS3)
  
  # EMOS with lognormal IGN minimization (i.e., ML estimation)
  EMOS4.fit <- fitMOS(train.eD, model = "lognormal", control = controlMOSlognormal(scoringRule = c("log")))
  EMOS4 <- quantileForecast(EMOS4.fit, test.eD, quantiles = Q)
  EMOS4 <- round(apply(EMOS4, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  colnames(EMOS4) <- paste0("EMOS4.", Q*100)
  fcst <- bind_cols(fcst, EMOS4)
  
  # EMOS with csg0 CRPS minimization
  EMOS5.fit <- fitMOS(train.eD, model = "gev0")
  EMOS5 <- quantileForecast(EMOS5.fit, test.eD, quantiles = Q)
  EMOS5 <- round(apply(EMOS5, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  colnames(EMOS5) <- paste0("EMOS5.", Q*100)
  fcst <- bind_cols(fcst, EMOS5)
  
  # quantile regression
  QR.fit <- quantreg::rq(observations~., data = train.eD[,-51], tau = Q, method = "br")
  QR <- predict(QR.fit, test.eD[,-51])
  QR <- round(apply(QR, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  QR <- t(apply(QR, 1, sort))
  colnames(QR) <- paste0("QR.", Q*100)
  fcst <- bind_cols(fcst, QR)
  
  # quantile regression with lasso
  # note that lasso QR runs into singular matrix problem, so we randomize the columns of the (already sorted) input matrix, for both train and test set
  set.seed(1234)
  rand.cols <- sample(1:50)
  QRL.fit <- quantreg::rq(observations~., data = bind_cols(train.eD[,rand.cols], observations = train.eD$observations), tau = Q, method = "lasso")
  #QRL.fit <- quantreg::rq(observations~., data = as_tibble(t(apply(train.eD[,-51], 1, jitter))), tau = Q, method = "lasso")
  QRL <- predict(QRL.fit, bind_cols(test.eD[,rand.cols], observations = test.eD$observations))
  QRL <- round(apply(QRL, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  QRL <- t(apply(QRL, 1, sort))
  colnames(QRL) <- paste0("QRL.", Q*100)
  fcst <- bind_cols(fcst, QRL)
  
  # quantile regression neural networks with 4 hidden neuron
  QRNN <- matrix(NA, nrow = nrow(test.eD), ncol = length(Q))
  for(tau in 1:length(Q))
  {
    set.seed(1234)
    QRNN.fit <- qrnn.fit(x = data.matrix(train.eD[,-c(51,52)]), y = matrix(train.eD$observations), n.hidden=4, n.trials=1, iter.max=500, tau = Q[tau], lower = 0)
    QRNN[,tau] <- qrnn.predict(data.matrix(test.eD[,-c(51,52)]), QRNN.fit)
  }
  QRNN <- round(QRNN, 1) # QRNN supports lower bound, only apply rounding here
  QRNN <- t(apply(QRNN, 1, sort))
  colnames(QRNN) <- paste0("QRNN4.", Q*100)
  fcst <- bind_cols(fcst, QRNN)
  
  # quantile regression neural networks with 8 hidden neuron
  QRNN <- matrix(NA, nrow = nrow(test.eD), ncol = length(Q))
  for(tau in 1:length(Q))
  {
    set.seed(1234)
    QRNN.fit <- qrnn.fit(x = data.matrix(train.eD[,-c(51,52)]), y = matrix(train.eD$observations), n.hidden=8, n.trials=1, iter.max=500, tau = Q[tau], lower = 0)
    QRNN[,tau] <- qrnn.predict(data.matrix(test.eD[,-c(51,52)]), QRNN.fit)
  }
  QRNN <- round(QRNN, 1) # QRNN supports lower bound, only apply rounding here
  QRNN <- t(apply(QRNN, 1, sort))
  colnames(QRNN) <- paste0("QRNN8.", Q*100)
  fcst <- bind_cols(fcst, QRNN)
  
  # quantile regression forest
  set.seed(1234)
  QRF.fit <-quantregForest::quantregForest(x = train.eD[,-c(51,52)], y = train.eD$observations)
  QRF<- predict(QRF.fit, newdata = test.eD[,-c(51,52)], what = Q)
  QRF <- round(apply(QRF, 2, function(x) ifelse(x<0, 0, x)), 1) # remove unphysical values and apply rounding
  QRF <- t(apply(QRF, 1, sort))
  colnames(QRF) <- paste0("QRF.", Q*100)
  fcst <- bind_cols(fcst, QRF)
  
  #####################################
  # write file
  #####################################
  file.name <- paste0(loc$stn[i],"_PP_", yr.te[1], "_", yr.te[length(yr.te)], ".txt")
  write.table(fcst, file = file.name, quote = FALSE, row.names = FALSE, sep = "\t")
  
  setTxtProgressBar(pb,i)
}
close(pb)

