

#' Write JAGS Model 4 to File
#'
#' @param model_filepath File path and name for JAGS model to be written
#'
#' @return Text file containing JAGS model is saved to \code{model_filepath}
#' @export
#'
write_model_4 <- function(model_filepath = "./models/model_4.txt"){
    if(!"models" %in% list.files("./"))
        dir.create("./models")
    cat("
        model{
            # multiple year, island-specific density, common growth rate, random observer-effect on sigma, no groups model
            # Prior distributions for model parameters
            mu_obs ~ dnorm(0, 0.001)
            tau_obs <- pow(sd_obs,-2)
            sd_obs ~ dunif(0,10)
            for(k in 1:(K-1)){
                alpha0[k] ~ dnorm(mu_obs,tau_obs) # scale parameter, observer-specific intercepts
            }
            alpha0[K] <- 0
            r ~ dunif(-5,5)
            for(s in 1:S){
                log_lambda0[s] ~ dunif(0,10)
            }
            # Individual level model: observations and process
            for(t in 1:T){
                for(s in 1:S){
                    for(j in 1:J){
                        for(i in 1:M){
                            z[i,j,s,t] ~ dbern(psi[j,s,t]) # Data augmentation variables
                            d[i,j,s,t] ~ dunif(0, B) # Distance is uniformly distributed
                            log(sigma[i,j,s,t]) <- alpha0[observer[i,j,s,t]]
                            mu[i,j,s,t] <- z[i,j,s,t]*exp(-d[i,j,s,t]*d[i,j,s,t]/(2*sigma[i,j,s,t]*sigma[i,j,s,t])) #p dep on dist class
                            # here using the half normal detection function (Buckland et al. 2001)
                            y[i,j,s,t] ~ dbern(mu[i,j,s,t])
                        }
                    }
                }
            }

            # Model for population size of groups
            for(s in 1:S){
                log_lambda_pop[s,1] <- log_lambda0[s]
                for(j in 1:J){
                    Nj[j,s,1] ~ dpois(lambda_offset[j,s,1])
                    log(lambda_offset[j,s,1]) <- log_lambda_pop[s,1] + log(transect_area[j,s,1])
                    # psi is a derived parameter
                    psi[j,s,1] <- lambda_offset[j,s,1]/M
                }
                for(t in 2:T){
                    log_lambda_pop[s,t] <- log_lambda_pop[s,t-1] + (r * yearstep[s,t-1])
                    for(j in 1:J){
                        Nj[j,s,t] ~ dpois(lambda_offset[j,s,t])
                        log(lambda_offset[j,s,t]) <- log_lambda_pop[s,t] + log(transect_area[j,s,t])
                        # psi is a derived parameter
                        psi[j,s,t] <- lambda_offset[j,s,t]/M
                    }
                }
            }
            # Derived quantities
            for(s in 1:S){
                for(t in 1:T){
                    lambda_pop[s,t] <- exp(log_lambda_pop[s,t])
                    EN_region[s,t] <- lambda_pop[s,t] * region_area[s]
                }
            }
        }
        ",
        fill = TRUE, file = model_filepath)
}
