
source("./R/summary.fn.R")
#design <- "r_25pct"
design <- "r_est"
years <- c(2018,2028)
fig_path <- paste0("./output/power_2surveys_",design)
filenames <- list.files(paste0(fig_path,"/mcmc", full.names = TRUE))
if(design == "r_est"){
    library(jagsUI)
    load("./output/out4.RData")
    r  <- out4$mean$r #r_est
}
if(design == "r_25pct"){
    r  <- log(0.75^(1/diff(years)))
}
nreps <- length(filenames)
cv <- bias <- cover_95 <- onetail_95 <- cover_90 <- onetail_90 <- rep(0,nreps)
print(system.time({
    for(rep in 1:nreps){
        cat("doing", rep, " of ", nreps, " reps\n\n")
        fname <- filenames[rep]
        jags_rep <- get(load(fname))
        rm(list = paste0("mcmc",substr(fname, regexpr("mcmc_",fname)[1] + 5, regexpr(".RData",fname)[1]-1)))
        r_post <- jags_rep$sims.list$r
        nsamples <- length(r_post)
        r_est <- summary.fn(r_post)
        onetail_95[rep] <- sum(r_post < r)/nsamples < 0.05 # one-tailed test for at least 25% decline at 0.95 prob
        onetail_90[rep] <- sum(r_post < r)/nsamples < 0.1 # one-tailed test for at least 25% decline at 0.90 prob
        cv[rep] <- as.numeric(abs(r_est["sd"]/r_est["mean"]))
        bias[rep] <- as.numeric(r_est["mean"] - r)
        lims <- r_est[c("2.5%","97.5%")]
        cover_95[rep] <- lims[1] < r & lims[2] > r
        lims <- r_est[c("5%","95%")]
        cover_90[rep] <- lims[1] < r & lims[2] > r
    }
}))
#
sq_err <- (r - bias)^2
sq_err_sum <- summary.fn(sq_err)
cv_sum <- summary.fn(cv)
bias_sum <- summary.fn(bias)
rmse <- sqrt(mean(sq_err))
(ci_cover_95 <- sum(cover_95)/nreps)
(onetail_power_95 <- 1 - (sum(onetail_95)/nreps))
(ci_cover_90 <- sum(cover_90)/nreps)
(onetail_power_90 <- 1 - (sum(onetail_90)/nreps))
#
save(sq_err, cv, bias, rmse, ci_cover_95, ci_cover_90, onetail_power_95, onetail_power_90, file = paste0(fig_path, "/summaries.RData"))
#
png(paste0(fig_path, "/cv_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
hist(cv, xlim = c(0,0.4), breaks = seq(0,0.4,0.005), xlab = "Coefficient of variation", main = "")
abline(v = cv_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png(paste0(fig_path, "/bias_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
hist(bias, xlim = c(-0.05,0.05), breaks = seq(-0.05, 0.05, 0.002), xlab = "Bias", main = "")
abline(v = bias_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
#
png(paste0(fig_path, "/sq_err_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
hist(sq_err, xlim = c(0,0.002), breaks = seq(0, 0.002, 0.0001), xlab = "Squared error", main = "")
abline(v = sq_err_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = c("red",rep("blue",2)))
dev.off()
