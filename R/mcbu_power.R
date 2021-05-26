

#' MCBU Model 4 Power Analysis - Data and MCMC Simulator
#'
#' @param log_lambda0 Numeric vector of length 2 containing starting island-specific MCBU densities in square km. Default set to preliminary estimates for 2018)
#' @param r Numeric value specifying exponential growth rate \code{r}. Default value is set to 25% decline over 10 years.
#' @param mu_obs Numeric value specifying the mean sigma value across observers on the log scale. Default set to preliminary estimates for 2003 and 2018.
#' @param sd_obs Numeric value specifying the standard deviation of sigma values across observers on the log scale. Default set to preliminary estimates for 2003 and 2018.
#' @param survey_years Numeric vector specifying in which years surveys will be conducted.
#' @param mcbu_data_file File path and name for MCBU distance data. Default is set to internal package \code{data} folder
#' @param save_output Logical value specifying whether to save JAGS output. Default set to FALSE.
#' @param output_filepath File path for saving JAGS output. Default is set to NULL which creates an internal \code{output} and related folders to store output.
#' @param batch Integer vector specifying replicates for data and MCMC simulation. Default set to single replicate.
#' @param M Integer vector specifying the augmented number of MCBU groups per square kilometer.
#' @param nb Number of MCMC samples to discard as burn-in samples.
#' @param nt Integer specifying the thinning rate for MCMC samples.
#' @param nc Integer specifying the number of individual MCMC chains to simulate.
#' @param use_parallel Logical value specifying whether to use parallel computing when simulating >1 MCMC chain. See help for \code{jags} function in \code{jagsUI} package for more details.
#'
#' @return No output objects are returned.
#' @export
#'
mcbu_power <- function(log_lambda0 = c(4.1257,3.9400), r = log(0.75^(1/10)), mu_obs = -2.9769, sd_obs = 0.4470, survey_years = c(2018,2028),
                       mcbu_data_file = "./data/distdata_mcbu.csv", save_output = FALSE, output_filepath = NULL,
                       batch = 1, M = c(200,100), ni = 200, nb = 100, nt = 1, nc = 1, use_parallel = TRUE){
    if(!"models" %in% list.files("./"))
        dir.create("./models")
    write_model_4()
    scenario <- paste0("mcbu_power_",length(survey_years),"surveys_r",round(r,2))
    if(save_output){
        if(is.null(output_filepath)){
            if(!"output" %in% list.files("./"))
                dir.create("./output")
            if(!scenario %in% list.files("./output"))
                dir.create(paste0("./output/", scenario))
            if(!"mcmc" %in% list.files(paste0("./output/", scenario)))
                output_filepath = paste0("./output/", scenario,"/mcmc")
            dir.create(output_filepath)
        }
    }
    B <- 0.1 # half-width of survey transects in km
    sites <- c("HALL","STMA")
    S <- length(sites)
    Js <- c(14,34)
    T <- length(survey_years)
    yearstep <- matrix(rep(diff(survey_years),S),S,T)
    K <- 8 + 1
    Jmat <- matrix(rep(Js,T), nrow = S, ncol = T, byrow = FALSE, dimnames = list(sites,paste0("y",survey_years)))
    params4 <- c("mu_obs", "sd_obs", "alpha0", "r", "lambda_pop", "EN_region")
    inits <- function(){list(mu_obs = -3, sd_obs = 1, r = 0, log_lambda0 = rep(4,S))}
    mcbu_data = read.csv(mcbu_data_file) # MCBU detections only, both survey years, data cleaned and ready for MCBU analysis
    exp_names = c("Region.Label","Area","Sample.Label","Effort","species","distance","size","obs","date","julian","year","bird_lat","bird_long","notes")
    if(!all(exp_names %in% colnames(mcbu_data))){
        stop(paste0("Column names missing from MCBU distance data. Expected column names are:\n\n",
                    paste0(exp_names, collapse = ", ")))
    }
    colnames(mcbu_data) <- c("region_id","region_area","survey_number","transect_length","species","distance","group_size","observer",
                             "date","julian","year","latitude","longitude","notes")
    indz <- unlist(mapply(function(i,n){rep(i,each = n)}, i = seq_len(nrow(mcbu_data)),
                          n = as.integer(mcbu_data$group_size), SIMPLIFY = TRUE))
    mcbu_data <- mcbu_data[indz,]
    mcbu_data$transect_id <- paste0(mcbu_data$region_id,"_",mcbu_data$survey_number)
    mcbu_data$survey_number <- as.character(mcbu_data$survey_number)
    mcbu_data$distancekm <- mcbu_data$distance / 1000
    observer_levels = data.frame(observer = c("MAL","REG","SMM","JAJ","RMR","ARD","MND","BWR","MDR","DRR","SLW"),
                                 level = 11:1)
    mcbu_data$observer <- as.factor(match(mcbu_data$observer, observer_levels$observer))
    data2003 <- subset(mcbu_data, mcbu_data$year==2003)
    data2018 <- subset(mcbu_data, mcbu_data$year==2018)
    split_data2003 <- lapply(split(data2003, data2003$region_id, drop = TRUE),function(x)split(x, x$transect_id))
    split_data2018 <- lapply(split(data2018, data2018$region_id, drop = TRUE),function(x)split(x, x$transect_id))
    datalist <- list("y2003" = split_data2003, "y2018" = split_data2018)
    region_area <- range(mcbu_data$region_area)
    transects <- data2018[match(unique(data2018$transect_id), data2018$transect_id),c("region_id","transect_length")]
    transects <- split(transects, transects$"region_id")
    # loop simulations
    for(rep in batch){
        cat("doing rep ", rep, " in ", min(batch), " to ", max(batch), "\n\n")
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
            mcmc_rep <- jags(jags_data, inits, params4, "./models/model_4.txt", n.thin = nt,
                             n.chains = nc, n.burnin = nb, n.iter = ni, parallel = use_parallel)
            # save output
            if(save_output){
                mcmc_name <- paste("mcmc",rep,sep="")
                file_name <- paste(output_filepath, "/mcmc_",rep,".RData",sep="")
                assign(mcmc_name, mcmc_rep)
                save(list=mcmc_name, file=file_name)
            }
        }))
    }
}
