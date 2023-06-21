#################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
# load the required packages
libs <- c("tidyverse", "patchwork", "lubridate", "SolarData", "xtable", "ggfan")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/dyang/Dropbox/Working papers/Combine/Data/Processed" # dir for processed data
stns <- c("bon", "dra", "fpk", "gwn", "psu", "sxf", "tbl")
methods <- c("EMOS1", "EMOS2", "EMOS3", "EMOS4", "EMOS5", "QR", "QRL", "QRNN4", "QRNN8", "QRF", "SQA", "CQRA")
yr.tr <- c(2019) # training years
yr.te <- c(2020) # verification years
Q <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99) # a set of quantiles to evaluate
m <- length(Q) # number of ensemble members
source("/Users/dyang/Dropbox/Working papers/Combine/Code/Functions.r")
plot.size = 7; line.size = 0.15; point.size = 0.5; legend.size = 0.4; text.size = 1.5;
# function to calculate pinball loss
pinball <- function(y, f, tau)
{
  pbs <- matrix(NA, nrow = nrow(f), ncol = ncol(f))
  for(j in 1:length(tau))
  {
    q <- f[,j]
    pbs[,j] <- ifelse(y >= q, tau[j] * (y-q), (1-tau[j]) * (q-y))
  }
  mean(colMeans(pbs)) # average n samples, then m quantiles
}
#################################################################################

# read metadata of SURFRAD stations
loc <- SURFRAD.loc
# reorder rows of "loc", to follow the order of the stations
loc <- loc[match(stns, loc$stn),]

setwd(dir0)

#################################################################################
# Evaluate quantile score
#################################################################################
# set empty matrix to hold forecasts errors
QS <- matrix(NA, nrow = length(methods), ncol = length(stns))
# get empty data.frame to hold new form of data to plot reliability
data.reli = data.band = data.cover = data.sharp <- NULL
pb = txtProgressBar(min = 0, max = length(stns), initial = 0, style = 3) 
for(i in 1:length(stns))
{
  # read data and separate in-sample and out-of-sample data
  data <- tibble(read.table(file.path(dir0, paste0(stns[i],"_qr_2019_2020_new_third.txt")), sep = "\t", header = TRUE)) %>%
    mutate(Time = ymd_hms(Time))
  data.tr <- data %>% filter(year(Time) %in% yr.tr)
  data.te <- data %>% filter(year(Time) %in% yr.te)
  
  #####################################
  # out-of-sample error computation
  #####################################
  fcst <- data.te
  obs <- data.te$y
  for(j in 1:length(methods))
  {
    assign(paste0("f.", methods[j]), data.matrix(fcst[,which(colnames(fcst) %in% paste0(methods[j], ".", Q*100))]))
    QS[j, i] <- pinball(obs, get(paste0("f.", methods[j])), Q)
  }
  
  #####################################
  # construct new data.frame format for reliability
  #####################################
  for(j in 1:length(methods))
  {
    qtau <- data.matrix(bind_cols(get(paste0("f.", methods[j]))))
    obs <- data.te$y
    n <- length(obs)
    P1 <- array(NA, length(Q))
    quant <- seq(0.01, 0.99, by = 0.02)
    P.quant <-  matrix(NA, nrow = length(Q), ncol = length(quant))
    for(k in 1:length(Q))
    {
      P1[k] <- length(which(obs < qtau[,k]))/n
      P.quant[k,] <- qbinom(quant, n, Q[k], lower.tail = TRUE)/n
    }
    # make data frame for reliability line
    tmp <- tibble(x = Q, y = P1) %>%
      mutate(method = methods[j], stn = toupper(stns[i]))
    data.reli <- bind_rows(data.reli, tmp)
    # make data frame for consistency band
    rownames(P.quant) <- as.character(Q)
    colnames(P.quant) <- quant
    P.quant <- reshape::melt(P.quant)
    tmp <- tibble(P.quant) %>%
      rename(x = X1, quant = X2) %>%
      mutate(method = methods[j], stn = toupper(stns[i]))
    data.band <- bind_rows(data.band, tmp)
  }
  
  #####################################
  # construct new data.frame format for reliability (coverage)
  #####################################
  for(j in 1:length(methods))
  {
    tmp <- bind_cols(data.te[,1:2], get(paste0("f.", methods[j]))) %>%
      gather(., key = "quantile", value = "value", -c(Time, y)) %>%
      mutate(quantile = as.numeric(str_split(quantile, "\\.", simplify = TRUE)[,2])/100) %>%
      rename(truth = y) %>%
      mutate(method = methods[j], stn = toupper(stns[i]))
    data.cover <- bind_rows(data.cover, tmp)
  }
  
  #####################################
  # construct new data.frame format for sharpness
  #####################################
  for(j in 1:length(methods))
  {
    tmp <- get(paste0("f.", methods[j]))
    # get interval wisth at three nominal coverage rates
    int.width.98 <- tmp[, match(paste0(methods[j], ".99"), colnames(tmp))] - tmp[, match(paste0(methods[j], ".1"), colnames(tmp))]
    int.width.90 <- tmp[, match(paste0(methods[j], ".95"), colnames(tmp))] - tmp[, match(paste0(methods[j], ".5"), colnames(tmp))]
    int.width.80 <- tmp[, match(paste0(methods[j], ".90"), colnames(tmp))] - tmp[, match(paste0(methods[j], ".10"), colnames(tmp))]
    tmp <- bind_cols(data.te[,1], a = int.width.80, b = int.width.90, c = int.width.98)
    names(tmp)[2:4] <- c("80", "90", "98")
    tmp <- tmp %>%
      tidyr::gather(., key = "interval", value = "width", -Time) %>%
      mutate(method = methods[j], stn = toupper(stns[i]))
    data.sharp <- bind_rows(data.sharp, tmp)
  }

  setTxtProgressBar(pb,i)
}
close(pb)

# print error table
# CRPS
error1 <- matrix(as.character(format(round(QS, 1), nsmall = 0)), ncol = ncol(QS))
for(i in 1:ncol(QS))
{
  bf.ind <- which.min(abs(QS[,i]))
  error1[bf.ind,i] <- paste0("\\textbf{", error1[bf.ind,i], "}", sep = "")
}
#error_print <- cbind(matrix(c("EMOS-N-CRPS", "EMOS-N-IGN", "EMOS-logN-CRPS", "EMOS-logN-IGN", "EMOS-gev0", methods[6:length(methods)])), error1)
error_print <- cbind(matrix(methods), error1)
print(xtable(error_print), include.rownames=FALSE, sanitize.text.function = function(x) sanitize.numbers(str = x, type = "latex"), booktabs = TRUE)

#################################################################################
# Reliability
#################################################################################
data.band$method <- factor(data.band$method, levels=methods)
levels(data.band$method) <- methods
data.reli$method <- factor(data.reli$method, levels=methods)
levels(data.reli$method) <- methods

# not used, because this is the same as Gneiting's plot
ggplot() + 
  geom_fan(data = data.band, aes(x=x,y=value, quantile=quant)) +
  facet_grid(stn~method) +
  geom_point(data = data.reli, aes(x = x, y = y), size = point.size*2) +
  geom_line(data = data.reli, aes(x = x, y = y), size = line.size*2) +
  #geom_abline(slope = 1, intercept = 0, size = line.size*2, linetype = "solid") +
  scale_x_continuous(name = expression(paste("Nominal probability, ", italic(tau))), limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_y_continuous(name = expression(paste("Observed proportion, ", bar(italic(z))[italic(tau)])), limits = c(0,1)) +
  scale_fill_gradient2(name = "Confidence", limits = c(0,1), midpoint=0,low="gray30",mid="grey50",high="gray80", breaks=c(0.25, 0.5, 0.75), labels=c("25%", "50%", "75%"))+
  theme_minimal() +
  theme(plot.margin = unit(c(0.3,0.2,0,0.1), "lines"), panel.spacing = unit(0.05, "lines"), plot.background = element_rect(fill = "transparent", color = NA), text = element_text(family = "Times", size = plot.size), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "lines"), size = plot.size), strip.text.y = element_text(margin = margin(0,0.05,0,0.05, "lines"), size = plot.size), axis.title = element_text(size = plot.size), axis.text = element_text(size = plot.size), legend.position = "bottom", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(0.5, "lines"), legend.key.width = unit(1.4, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", colour = "transparent")) # 


#################################################################################
# Reliability (coverage)
#################################################################################
data.cover2 <- data.cover %>%
  filter(week(Time) == 1)

# consistency bands
coverage.plot <- data.cover %>%
  group_by(stn, method, quantile) %>%
  coverage(., date_column = Time, band_type = "consistency")
coverage.plot$method <- factor(coverage.plot$method, levels=methods)
levels(coverage.plot$method) <- methods

p1 <- ggplot(coverage.plot) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.3, linetype = "solid", colour = "grey70") +
  geom_errorbar(aes(x = quantile, ymin = l, ymax = u), width = 0.05, size = 0.3, colour = "black") +
  geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.7) +
  geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.3) +
  #geom_line(aes(x = quantile, y = l), size = line.size*2)+
  facet_grid(stn~method) +
  scale_x_continuous(name = "Quantile level", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))+
  scale_y_continuous(name = "Coverage", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.05, "lines"), plot.background = element_rect(fill = "transparent", color = NA), text = element_text(family = "Times", size = plot.size), strip.text.x = element_text(margin = ggplot2::margin(0.05,0,0.05,0, "lines"), size = plot.size), strip.text.y = element_text(margin = ggplot2::margin(0,0.05,0,0.05, "lines"), size = plot.size), axis.title = element_text(size = plot.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(1, "lines"), legend.key.width = unit(0.9, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.key = element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent", colour = "transparent"), panel.grid.minor = element_blank())

setwd("/Users/dyang/Dropbox/Working papers/Combine/tex")
ggsave("reli.pdf", p1, device = cairo_pdf, width = 180, height = 120, unit = "mm")

#################################################################################
# Sharpness
#################################################################################
data.sharp$method <- factor(data.sharp$method, levels=methods)

p2 <- ggplot(data.sharp) +
  #geom_violin(aes(y = width, x = interval, fill = method), size = line.size, scale = "width") +
  geom_boxplot(aes(y = width, x = interval, fill = interval), size = line.size, outlier.size = point.size, outlier.stroke = 0, outlier.alpha = 0.1, outlier.shape = '.', alpha = 0.7) +
  facet_grid(stn~method) + 
  scale_fill_manual(name = "", values = ggthemes::colorblind_pal()(8)[2:4]) +
  scale_x_discrete(name = "Nominal coverage rate [%]") +
  scale_y_continuous(name = expression(paste("Interval width [W/", m^2, "]")), limits = c(0,1000), breaks = c(0, 400, 800)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.05, "lines"), plot.background = element_rect(fill = "transparent", color = NA), text = element_text(family = "Times New Roman", size = plot.size), strip.text.x = element_text(margin = ggplot2::margin(0.05,0,0.05,0, "lines"), size = plot.size), strip.text.y = element_text(margin = ggplot2::margin(0,0.05,0,0.05, "lines"), size = plot.size), axis.title = element_text(size = plot.size), axis.text = element_text(size = plot.size), legend.position = "none", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(1, "lines"), legend.key.width = unit(0.9, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), panel.background = element_rect(fill = "transparent", colour = "transparent")) 

setwd("/Users/dyang/Dropbox/Working papers/Combine/tex")
ggsave("sharp.pdf", p2, device = cairo_pdf, width = 180, height = 120, unit = "mm")


#################################################################################
# table for model weights
#################################################################################
weights <- tibble(read.table(file.path(dir0, paste0(stns[4],"_qr_2019_2020_Wnq.txt")), sep = "\t", header = TRUE)) 
table_print <- cbind(matrix(methods[1:10]), weights)
print(xtable(table_print), include.rownames=FALSE, sanitize.text.function = function(x) sanitize.numbers(str = x, type = "latex"), booktabs = TRUE)







# #####################################
# # make reference forecast
# #####################################
# # simple quantile averaging
# f.SQA <- matrix(NA, nrow = nrow(fcst), ncol = length(Q))
# for(j in 1:length(Q))
# {
#   f.SQA[,j] <- rowMeans(cbind(f.QR[,j], f.QRL[,j], f.QRNN4[,j], f.QRNN8[,j], f.QRF[,j]))
# }
# colnames(f.SQA) <- paste0("SQA.", Q*100)
# QS.out[6, i] <- pinball(obs, f.SQA, Q)


