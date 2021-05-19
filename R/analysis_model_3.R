# load packages and functions
library(jagsUI)
source("./R/summary.fn.R")

# process data
distdata_mcbu <- read.csv("./data/distdata_mcbu.csv") # MCBU detections only, both survey years, data cleaned and ready for MCBU analysis
colnames(distdata_mcbu) <- c("region_id","region_area","survey_number","transect_length","species","distance","group_size","observer",
                             "date","julian","year","latitude","longitude","notes")
distdata_mcbu$transect_id <- paste0(distdata_mcbu$region_id,"_",distdata_mcbu$survey_number)
distdata_mcbu$survey_number <- as.character(distdata_mcbu$survey_number)
distdata_mcbu$distancekm <- distdata_mcbu$distance / 1000
# Rename observer levels in the dataframe to reflect ordinal rank:
observer_levels <- read.csv("./data/observer_factor_levels.csv")
distdata_mcbu$observer <- as.factor(match(distdata_mcbu$observer, observer_levels$observer))
data2003 <- subset(distdata_mcbu, distdata_mcbu$year==2003)
data2018 <- subset(distdata_mcbu, distdata_mcbu$year==2018)
split_data2003 <- lapply(split(data2003, data2003$region_id, drop = TRUE),function(x)split(x, x$transect_id))
split_data2018 <- lapply(split(data2018, data2018$region_id, drop = TRUE),function(x)split(x, x$transect_id))
datalist <- list("y2003" = split_data2003, "y2018" = split_data2018)
#
yearstep <- matrix(c(15,15),2,1)
region_area <- range(distdata_mcbu$region_area)
K <- nrow(observer_levels) + 1
B <- 0.1 # half-width of survey transects in km
M <- c(200,100)
J <- sapply(datalist, function(x)sapply(x,length))
S <- nrow(J)
T <- ncol(J)
y <- z <- group_size <- d <- observer <- array(0,c(max(M),max(J),S,T))
z <- group_size <- d <- observer <- array(NA,c(max(M),max(J),S,T))
transect_area <- array(0,c(max(J),S,T))
for(t in 1:T){
    for(s in 1:S){
        for(j in 1:J[s,t]){
            nobs <- nrow(datalist[[t]][[s]][[j]])
            y[1:nobs,j,s,t] <- z[1:nobs,j,s,t] <- 1
            group_size[1:nobs,j,s,t] <- datalist[[t]][[s]][[j]]$group_size - 1 # subtract one to use zero-truncated Poisson distribution for group size
            d[1:nobs,j,s,t] <- datalist[[t]][[s]][[j]]$distancekm
            transect_area[j,s,t] <- datalist[[t]][[s]][[j]]$transect_length[1] * B * 2
            observer[,j,s,t] <- as.integer(datalist[[t]][[s]][[j]]$observer[1])
        }
    }
    z[,(min(J[,t])+1):max(J[,t]),1,t] <- 0
    z[(min(M)+1):max(M),,1,t] <- 0
    observer[,(min(J[,t])+1):max(J[,t]),1,t] <- K
}
# format jags input data
jags.data <- list(M = max(M),
                  J = max(J),
                  S = S,
                  T = T,
                  K = K,
                  y = y,
                  z = z,
                  d = d,
                  group_size = group_size,
                  observer = observer,
                  transect_area = transect_area,
                  region_area = region_area,
                  B = B,
                  yearstep = yearstep)
# set up, run jags, and save output
ni <- 105000
nb <- 5000
nt <- 3
nc <- 3
params3 <- c("mu_obs", "sd_obs", "alpha0", "alpha1", "r", "lambda_pop", "lambda_group", "lambda_grsize", "EN_region")
inits <- function(){list(mu_obs = -3, sd_obs = 1, alpha1 = 0, r = 0, log_lambda0 = rep(4,S), lambda_group = rep(1,S))}
out3 <- jags(jags.data, inits, params3, "multi_year_open_pop_observerRE_model.txt", n.thin = nt,
             n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
save(out3full = out3, file = "./output/out3.RData")


traceplot(out3)
options(max.print = 50000)
print(out3, 3)


# derive abundances
EN_total <- t(apply(out3$sims.list$EN_region,1,colSums))

# make some figures
png("./output/EN_total_model3.png",width = 6.5, height = 9, units = "in", res = 192)
par(mfrow=c(2,1))
xs <- summary.fn(EN_total[,1])[c("mean","2.5%","97.5%")]
hist(EN_total[,1], xlim = c(10000,40000), breaks = seq(10000,40000,500),
     main = paste0("Total abundance 2003 (N-hat = ",round(xs[1],0),")"), xlab = "Abundance")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
xs <- summary.fn(EN_total[,2])[c("mean","2.5%","97.5%")]
hist(EN_total[,2], xlim = c(10000,40000), breaks = seq(10000,40000,500),
     main = paste0("Total abundance 2003 (N-hat = ",round(xs[1],0),")"), xlab = "Abundance")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/15yr_pctdecline_total_model3.png",width = 6.5, height = 9, units = "in", res = 192)
pctdecline <- 100 * (1 - (EN_total[,2]/EN_total[,1]))
xs <- summary.fn(pctdecline)[c("mean","2.5%","97.5%")]
hist(pctdecline, xlim = c(20,60), breaks = seq(20,60,1),
     main = paste0("15-year percent decline (pct-hat = ",round(xs[1],1),")"), xlab = "Percent decline")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/growth_rate_r_model3.png",width = 6.5, height = 6.5, units = "in", res = 192)
xs <- out3$summary[rownames(out3$summary)=="r",colnames(out3$summary)%in%c("mean","2.5%","97.5%")]
hist(out3$sims.list$r, xlim = c(-0.05,-0.015), breaks = seq(-0.05,-0.015,0.001),
     main = paste0("Exponential growth rate (r-hat = ",round(xs[1],4),")"), xlab = "Growth rate")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/growth_rate_lambda_model3.png",width = 6.5, height = 6.5, units = "in", res = 192)
grlambda <- exp(out3$sims.list$r)
xs <- summary.fn(grlambda)[c("mean","2.5%","97.5%")]
hist(grlambda, xlim = c(0.90,1.05), breaks = seq(0.90,1.05,0.001),
     main = paste0("Geometric growth rate (lambda-hat = ",round(xs[1],4),")"), xlab = "Growth rate")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/10yr_prjN_total_model3.png",width = 9, height = 6.5, units = "in", res = 192)
prjN <- EN_total[,2] * exp(out3$sims.list$r)^10
xs1 <- summary.fn(prjN)[c("mean","2.5%","97.5%")]
hist(prjN, xlim = c(5000,40000), breaks = seq(5000,40000,500), col = rgb(211/255,211/255,211/255,0.75), ylim = c(0,30000),
     main = paste0("Projected total abundance 2028 (N-hat = ",round(xs1[1],0),")"), xlab = "Abundance")
xs2 <- summary.fn(EN_total[,2])[c("mean","2.5%","97.5%")]
hist(EN_total[,2], breaks = seq(5000,40000,500), col = rgb(211/255,211/255,211/255,0.75), add = T)
xs3 <- summary.fn(EN_total[,1])[c("mean","2.5%","97.5%")]
hist(EN_total[,1], breaks = seq(5000,40000,500), col = rgb(211/255,211/255,211/255,0.75), add = T)
abline(v = c(xs1[2:3],xs2[2:3],xs3[2:3]), lty = 2, lwd = 2, col = "blue")
segments(x0 = c(xs1[1],xs2[1],xs3[1]), y0 = rep(0,3), x1 = c(xs1[1],xs2[1],xs3[1]), y1 = rep(28000,3),
         lty = 1, lwd = 2, col = "red")
text(c(xs1[1],xs2[1],xs3[1]), 30000, c("2028","2018","2003"), cex=0.75)
dev.off()
