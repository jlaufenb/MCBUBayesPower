# run this script by using the "run_power.R" script
library(jagsUI)
source("./R/summary.fn.R")
B <- 0.1 # half-width of survey transects in km
sites <- c("HALL","STMA")
S <- length(sites)
M <- c(200,100)
Js <- c(14,34)
years <- c(2018,2028)
T <- length(years)
yearstep <- matrix(rep(diff(years),S),S,T)
K <- 8 + 1
Jmat <- matrix(rep(Js,2), nrow = S, ncol = T, byrow = FALSE, dimnames = list(sites,paste0("y",years)))
ni <- 5500
nb <- 500
nt <- 1
nc <- 1
params4 <- c("mu_obs", "sd_obs", "alpha0", "r", "lambda_pop", "EN_region")
# Call JAGS (ART 1.4 min), check convergence and summarize posterior distributions
inits <- function(){list(mu_obs = -3, sd_obs = 1, r = 0, log_lambda0 = rep(4,S))}
if(!"power_2surveys_r_est" %in% list.files("./output"))
    dir.create("./output/power_2surveys_r_est")
if(!"mcmc" %in% list.files("./output/power_2surveys_r_est"))
    dir.create("./output/power_2surveys_r_est/mcmc")
load("./output/out4.RData")
log_lambda0 <- c(summary.fn(log(out4$sims.list$lambda_pop[,1,2]))["mean"],
                 summary.fn(log(out4$sims.list$lambda_pop[,2,2]))["mean"]) #log_lambda0_est #vector of length S
r  <- out4$mean$r #r_est
mu_obs <- out4$mean$mu_obs #mu_obs_est
sd_obs <- out4$mean$sd_obs #sd_obs_est
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
region_area <- range(distdata_mcbu$region_area)
transects <- data2018[match(unique(data2018$transect_id), data2018$transect_id),c("region_id","transect_length")]
transects <- split(transects, transects$"region_id")
# loop simulations
for(rep in batch){
    cat("doing", rep, " of ", length(batch), " reps\n\n")
    print(system.time({
        set.seed(2003 + rep)
        #
        alpha <- rep(0,K-1)
        log_lambda_pop <- lambda_pop <- EN_region <- matrix(0,S,T)
        observer <- array(K,c(max(M),max(Jmat),S,T))
        y <- z <- d <- sigma <- mu <- array(0,c(max(M),max(Jmat),S,T))
        transect_area <- array(0,c(max(Jmat),S,T))
        lambda_offset <- array(0,c(max(Jmat),S,T))
        psi <- array(0,c(max(Jmat),S,T))
        alpha0 <- c(rnorm(K-1,mu_obs,sd_obs),0) # scale parameter, observer-specific intercepts
        # Individual level model: observations and process
        for(s in 1:S){
            log_lambda_pop[s,1] <- log_lambda0[s]
            for(j in 1:Jmat[s,1]){
                transect_area[j,s,1] <- transects[[s]][j,"transect_length"] * B * 2
                lambda_offset[j,s,1] <- exp(log_lambda_pop[s,1]) * transect_area[j,s,1]
                # psi is a derived parameter
                psi[j,s,1] <- lambda_offset[j,s,1]/M[1]
            }
            for(t in 2:T){
                log_lambda_pop[s,t] <- log_lambda_pop[s,t-1] + (r * yearstep[s,t-1])
                for(j in 1:Jmat[s,t]){
                    transect_area[j,s,t] <- transects[[s]][j,"transect_length"] * B * 2
                    lambda_offset[j,s,t] <- exp(log_lambda_pop[s,t]) * transect_area[j,s,t]
                    # psi is a derived parameter
                    psi[j,s,t] <- lambda_offset[j,s,t]/M[t]
                }
            }
            for(t in 1:T){
                lambda_pop[s,t] <- exp(log_lambda_pop[s,t])
                EN_region[s,t] <- lambda_pop[s,t] * region_area[s]
                for(j in 1:Jmat[s,t]){
                    obsvr <- sample(1:(K-1),1)
                    observer[,j,s,t] <- obsvr
                    zs <- rbinom(1,M[t],psi[j,s,t])
                    z[1:zs,j,s,t] <- 1
                    d[,j,s,t] <- runif(M[t],0,B) # Distance is uniformly distributed
                    mu[,j,s,t] <- z[,j,s,t]*exp(-d[,j,s,t]*d[,j,s,t]/(2*exp(alpha0[obsvr])*exp(alpha0[obsvr]))) #p dep on dist class
                    y[1:M[t],j,s,t] <- rbinom(M[t],1,mu[,j,s,t])
                    z[y[,j,s,t]==0,j,s,t] <- NA
                    d[y[,j,s,t]==0,j,s,t] <- NA
                }
                if(Jmat[s,t]<max(Jmat)){
                    z[,(Jmat[s,t]+1):max(Jmat),s,t] <- 0
                }
                if(M[t]<max(M)){
                    z[(min(M[t])+1):max(M[t]),,s,t] <- 0
                }
            }
        }
        # Compile jags data list
        jags_data <- list(M = max(M),
                          J = max(Jmat),
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
        # run mcmc
        mcmc_rep <- jags(jags_data, inits, params4, "multi_year_open_pop_observerRE_nogroup_model.txt", n.thin = nt,
                         n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
        # save output
        mcmc_name <- paste("mcmc",rep,sep="")
        file_name <- paste("./output/power_2surveys_r_est/mcmc/mcmc_",rep,".RData",sep="")
        assign(mcmc_name, mcmc_rep)
        save(list=mcmc_name, file=file_name)
    }))
}
