#################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
# load the required packages
libs <- c("dplyr", "lubridate", "SolarData", "ensembleMOS", "xtable")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/dyang/Dropbox/Working papers/Combine/Data/Processed" # dir for processed data
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
# objective function for IGN minimization
objectiveFUN <- function(pars) {
  if (control$coefRule == "square") {
    pars[-c(1:3)] <- pars[-c(1:3)]^2
  }
  if (control$varRule == "square") {
    pars[1] <- pars[1]^2
  }
  mval <- xTrain %*% pars[-c(1:2)]
  vval <- pars[1] + pars[2]^2 * var
  mu <- log(mval^2/sqrt(vval + mval^2))
  sig.sq <- log(1 + vval/mval^2)
  ign <- log(obs) + 0.5 * log(2 * pi * sig.sq) + ((log(obs) - mu)^2)/(2 * sig.sq)
  ign <- ign[!(obs == 0)]
  sum(ign)
}
#################################################################################

# read metadata of SURFRAD stations
loc <- SURFRAD.loc
# reorder rows of "loc", to follow the order of the stations
loc <- loc[match(stns, loc$stn),]

#################################################################################
# Check GWN EMOS4 forecasts
#################################################################################
setwd(dir0)

data <- tibble(read.table(file.path(dir0, paste0(stns[4],"_ECMWF_2017_2020.txt")), sep = "\t", header = TRUE)) %>%
  mutate(Time = ymd_hms(Time)) %>%
  mutate(Time = Time + loc$tz[4]*3600) %>%
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


# EMOS with lognormal IGN minimization (i.e., ML estimation)
EMOS4.fit <- fitMOS(train.eD, model = "lognormal", control = controlMOSlognormal(scoringRule = c("log"), optimRule = "BFGS")) # Nelder-Mead or BFGS
# fitting v, m, sigma, mu
v <- EMOS4.fit$c + EMOS4.fit$d*apply(train.eD[,1:50], 1, var)
m <- EMOS4.fit$a + data.matrix(train.eD[,1:50])%*%matrix(EMOS4.fit$B)
sigma <- sqrt(log(1+v/m^2))
mu <- log(m^2/sqrt(v+m^2))
# IGN
obj <- sum(log(train.eD$observations) + 0.5 * log(2 * pi * sigma^2) + ((log(train.eD$observations) - mu)^2)/(2 * sigma^2))
pars <- tibble(beta1 = EMOS4.fit$a, beta2 = EMOS4.fit$B[1], a = EMOS4.fit$c, b = EMOS4.fit$d, obj = obj)

# EMOS with lognormal IGN minimization (i.e., ML estimation)
EMOS4.fit <- fitMOS(train.eD, model = "lognormal", control = controlMOSlognormal(scoringRule = c("log"), optimRule = "Nelder-Mead")) # Nelder-Mead or BFGS
# fitting v, m, sigma, mu
v <- EMOS4.fit$c + EMOS4.fit$d*apply(train.eD[,1:50], 1, var)
m <- EMOS4.fit$a + data.matrix(train.eD[,1:50])%*%matrix(EMOS4.fit$B)
sigma <- sqrt(log(1+v/m^2))
mu <- log(m^2/sqrt(v+m^2))
# IGN
obj <- sum(log(train.eD$observations) + 0.5 * log(2 * pi * sigma^2) + ((log(train.eD$observations) - mu)^2)/(2 * sigma^2))
pars <- bind_rows(pars, tibble(beta1 = EMOS4.fit$a, beta2 = EMOS4.fit$B[1], a = EMOS4.fit$c, b = EMOS4.fit$d, obj = obj))
pars


#######################
# no longer used
#######################

EMOS4 <- quantileForecast(EMOS4.fit, test.eD, quantiles = Q)

# second term
0.5 * log(2 * pi * sigma^2)

# out-of-sample prediction
v <- EMOS4.fit$c + EMOS4.fit$d*apply(test.eD[,1:50], 1, var)
m <- EMOS4.fit$a + data.matrix(test.eD[,1:50])%*%matrix(EMOS4.fit$B)
sigma <- sqrt(log(1+v/m^2))
mu <- log(m^2/sqrt(v+m^2))
qlnorm(Q, meanlog = mu[3], sdlog = sigma[3])


threshold <- seq(0,6000, by=1)
predDist <- matrix(NA, nrow = length(threshold), ncol = 10)
for(k in 1:10)
{
  predDist[,k] <- dlnorm(threshold, meanlog = mu[k], sdlog = sigma[k])
}
plot(threshold, predDist[,3], type ="l")




# trace the optimization
ensembleData <-train.eD
control <- controlMOSlognormal(scoringRule = c("log"))

na.rm <- FALSE
forecasts <- ensembleForecasts(ensembleData)
exchangeable <- ensembleGroups(ensembleData)

"varNA" <- function(x) {
  v <- var(x, na.rm = TRUE)
}
nEX <- 2
if (!(nullX <- is.null(exchangeable))) {
  "catStrings" <- function(strings) {
    if (length(strings) == 1) 
      return(strings)
    if (is.null(strings[1]) || !is.character(strings)) {
      stop("Not a valid string")
    }
    newstrings <- strings[1]
    for (i in 2:length(strings)) {
      newstrings <- paste(newstrings, strings[i], sep = ".")
    }
    newstrings
  }
  namEX <- as.character(exchangeable)
  uniqueEX <- unique(namEX)
  nEX <- length(uniqueEX)
  splitEX <- split(seq(along = exchangeable), exchangeable)
  for (i in 1:nEX) {
    if (length(splitEX[[i]]) > 1) {
      M <- apply(forecasts[, splitEX[[i]]], 1, function(z) all(is.na(z)))
    }
  }
  D <- apply(forecasts, 1, function(z) !any(is.na(z)))
  nObs <- dim(forecasts[D, ])[1]
  matEX <- matrix(NA, nObs, ncol = nEX)
  for (k in 1:nEX) {
    if (length(splitEX[[k]]) > 1) {
      matEX[, k] <- apply(forecasts[D, splitEX[[k]]], 
                          1, mean)
    }
    else {
      matEX[, k] <- forecasts[D, splitEX[[k]]]
    }
  }
  ensembleEX <- cbind(matEX, dates = ensembleData$dates[D], 
                      observations = ensembleData$obs[D])
}
fcstEX <- data.frame(matEX)
fcstHour <- ensembleFhour(ensembleData)
initTime <- ensembleItime(ensembleData)
fitData <- ensembleData(forecasts = fcstEX, dates = ensembleData$dates[D], 
                        observations = ensembleData$obs[D], forecastHour = fcstHour, 
                        initializationTime = initTime)
fitForecasts <- ensembleForecasts(fitData)
obs <- fitData$observations
xTrain <- cbind(rep(1, length(obs)), fitForecasts)
var <- apply(ensembleForecasts(ensembleData)[D, ], 1, varNA)
if (any(is.null(c(control$start$c, control$start$d)))) {
  control$start$c <- 5
  control$start$d <- 1
}
if (any(is.null(c(control$start$B, control$start$a)))) {
  olsCoefs <- lm(obs ~ ensembleForecasts(fitData))$coef
  control$start$a <- olsCoefs[1]
  B <- olsCoefs[-1]
  if (control$coefRule == "square") {
    control$start$B <- sqrt(abs(B))
  }
  else control$start$B <- B
}

c <- control$start$c
if (control$varRule == "square") 
  c <- sqrt(c)
pars <- c(c, sqrt(control$start$d), control$start$a, control$start$B)

control$optimRule <- "Nelder-Mead"
trace(objectiveFUN, exit = quote(print(as.numeric(pars))))
out <- capture.output(opt <- optim(pars, fn = objectiveFUN, method = control$optimRule, 
             control = list(maxit = control$maxIter)))

optB <- opt$par[-c(1:3)]
B <- rep(NA, nForecasts)
for (i in 1:nEX) {
  B[splitEX[[i]]] <- optB[i]
}
if (control$coefRule == "square") {
  B <- B^2
  if (!nullX) {
    for (i in 1:nEX) {
      B[splitEX[[i]]] <- B[splitEX[[i]]]/length(splitEX[[i]])
    }
  }
  if (control$varRule == "square") {
    opt$par[1] <- opt$par[1]^2
  }
  fit <- structure(list(a = opt$par[3], B = B, c = opt$par[1], 
                        d = opt$par[2]^2, exhangeable = exchangeable), class = "fitMOSlognormal")
}

fit

# pick out the parameters
tmp <- lapply(strsplit(out[(1:(length(out)/2))*2], split = " +"), function(x) x[c(2:5)])
tmp <- t(sapply(tmp, as.numeric))
colnames(tmp) <- c("c", "d", "a", "B")
tmp <- as_tibble(tmp) %>%
  mutate(B = B/50) %>%
  mutate(d = d^2)
# calculate IGN in each iteration
IGN <- array(NA, nrow(tmp))
S2 <- apply(train.eD[,1:50], 1, var)
Xsum <- rowSums(data.matrix(train.eD[,1:50]))
for(i in 1:nrow(tmp))
{
  par <- tmp[i,]
  v <- par$c + par$d*S2
  m <- par$a + Xsum*par$B
  sigma <- sqrt(log(1+v/m^2))
  mu <- log(m^2/sqrt(v+m^2))
  # IGN
  ign <- log(train.eD$observations) + 0.5 * log(2 * pi * sigma^2) + ((log(train.eD$observations) - mu)^2)/(2 * sigma^2)
  IGN[i] <- sum(ign)
}

plot(tmp$c, ylim = c(0,1.6e+7))
plot(tmp$d, ylim = c(0,1.6e+7))

