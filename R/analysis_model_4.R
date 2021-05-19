# load packages and functions
library(jagsUI)
source("./R/summary.fn.R")

# process data
distdata_mcbu <- read.csv("./data/distdata_mcbu.csv") # MCBU detections only, both survey years, data cleaned and ready for MCBU analysis
colnames(distdata_mcbu) <- c("region_id","region_area","survey_number","transect_length","species","distance","group_size","observer",
                             "date","julian","year","latitude","longitude","notes")
indz <- unlist(mapply(function(i,n){rep(i,each = n)}, i = seq_len(nrow(distdata_mcbu)),
                      n = as.integer(distdata_mcbu$group_size), SIMPLIFY = TRUE))
distdata_mcbu <- distdata_mcbu[indz,]
hist(distdata_mcbu$distance)
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
y <- z <- d <- observer <- array(0,c(max(M),max(J),S,T))
z <- d <- observer <- array(NA,c(max(M),max(J),S,T))
transect_area <- array(0,c(max(J),S,T))
for(t in 1:T){
    for(s in 1:S){
        for(j in 1:J[s,t]){
            nobs <- nrow(datalist[[t]][[s]][[j]])
            y[1:nobs,j,s,t] <- z[1:nobs,j,s,t] <- 1
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
params4 <- c("mu_obs", "sd_obs", "alpha0", "r", "lambda_pop", "EN_region")
inits <- function(){list(mu_obs = -3, sd_obs = 1, r = 0, log_lambda0 = rep(4,S))}
out4 <- jags(jags.data, inits, params4, "multi_year_open_pop_observerRE_nogroup_model.txt", n.thin = nt,
             n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
save(out4full = out4, file = "./output/out4.RData")


traceplot(out4)
options(max.print = 50000)
print(out4, 3)


# derive abundances
EN_regions <- out4$sims.list$EN_region
EN_total <- t(apply(EN_regions,1,colSums))

# make some figures
png("./output/EN_regions_model4.png",width = 6.5, height = 9, units = "in", res = 192)
par(mfrow = c(2,1))
xs1 <- summary.fn(EN_regions[,1,1])[c("mean","2.5%","97.5%")]
hist(EN_regions[,1,1], xlim = c(500,2500), breaks = seq(500,2500,25), col = rgb(211/255,211/255,211/255,0.75),
     ylim = c(0,25000), main = "Halls Island abundance", xlab = "Abundance")
#
xs2 <- summary.fn(EN_regions[,1,2])[c("mean","2.5%","97.5%")]
hist(EN_regions[,1,2], breaks = seq(500,2500,25), col = rgb(211/255,211/255,211/255,0.75), add = T)
abline(v = c(xs1[2:3],xs2[2:3]), lty = 2, lwd = 2, col = "blue")
segments(x0 = c(xs1[1],xs2[1]), y0 = rep(0,2), x1 = c(xs1[1],xs2[1]), y1 = rep(28000,2),
         lty = 1, lwd = 2, col = "red")
text(c(xs1[1],xs2[1])-400, 20000, paste0(c("2003 N-hat = ","2018 N-hat = "), round(c(xs1[1],xs2[1]),0)), cex=0.75)
#
xs3 <- summary.fn(EN_regions[,2,1])[c("mean","2.5%","97.5%")]
hist(EN_regions[,2,1], xlim = c(10000,30000), breaks = seq(10000,30000,250), col = rgb(211/255,211/255,211/255,0.75),
     ylim = c(0,25000), main = "St Matthews Island abundance", xlab = "Abundance")
#
xs4 <- summary.fn(EN_regions[,2,2])[c("mean","2.5%","97.5%")]
hist(EN_regions[,2,2], breaks = seq(10000,30000,250), col = rgb(211/255,211/255,211/255,0.75), add = T)
abline(v = c(xs3[2:3],xs4[2:3]), lty = 2, lwd = 2, col = "blue")
segments(x0 = c(xs3[1],xs4[1]), y0 = rep(0,2), x1 = c(xs3[1],xs4[1]), y1 = rep(28000,2),
         lty = 1, lwd = 2, col = "red")
text(c(xs3[1],xs4[1])-4000, 20000, paste0(c("2003 N-hat = ","2018 N-hat = "), round(c(xs3[1],xs4[1]),0)), cex=0.75)
dev.off()
#
png("./output/EN_total_model4.png",width = 6.5, height = 9, units = "in", res = 192)
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
png("./output/15yr_pctdecline_total_model4.png",width = 6.5, height = 9, units = "in", res = 192)
pctdecline <- 100 * (1-(EN_total[,2]/EN_total[,1]))
xs <- summary.fn(pctdecline)[c("mean","2.5%","97.5%")]
hist(pctdecline, xlim = c(20,60), breaks = seq(20,60,1),
     main = paste0("15-year percent decline (pct-hat = ",round(xs[1],1),")"), xlab = "Percent decline")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/growth_rate_r_model4.png",width = 6.5, height = 6.5, units = "in", res = 192)
xs <- out4$summary[rownames(out4$summary)=="r",colnames(out4$summary)%in%c("mean","2.5%","97.5%")]
hist(out4$sims.list$r, xlim = c(-0.05,-0.015), breaks = seq(-0.05,-0.015,0.001),
     main = paste0("Exponential growth rate (r-hat = ",round(xs[1],4),")"), xlab = "Growth rate")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/growth_rate_lambda_model4.png",width = 6.5, height = 6.5, units = "in", res = 192)
grlambda <- exp(out4$sims.list$r)
xs <- summary.fn(grlambda)[c("mean","2.5%","97.5%")]
hist(grlambda, xlim = c(0.90,1.05), breaks = seq(0.90,1.05,0.001),
     main = paste0("Geometric growth rate (lambda-hat = ",round(xs[1],4),")"), xlab = "Growth rate")
abline(v = xs, lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png("./output/10yr_prjN_total_model4.png",width = 9, height = 6.5, units = "in", res = 192)
prjN <- EN_total[,2] * exp(out4$sims.list$r)^10
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
