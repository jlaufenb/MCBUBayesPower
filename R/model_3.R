

#' Write JAGS Model 3 to File
#'
#' @param model_filepath File path and name for JAGS model to be written
#'
#' @return Text file containing JAGS model is saved to \code{model_filepath}
#' @export
#'
write_model_3 <- function(model_filepath = "./models/model_3.txt"){
    cat("
        model{
            # multiple year, island-specific density and groups, common growth rate, random observer-effect on sigma model
            # Prior distributions for model parameters
            mu_obs ~ dnorm(0, 0.001)
            tau_obs <- pow(sd_obs,-2)
            sd_obs ~ dunif(0,10)
            for(k in 1:(K-1)){
                alpha0[k] ~ dnorm(mu_obs,tau_obs) # scale parameter, observer-specific intercepts
            }
            alpha0[K] <- 0
            alpha1 ~ dunif(-5,5) # scale parameter, group-size effect
            r ~ dunif(-5,5)
            for(s in 1:S){
                log_lambda0[s] ~ dunif(0,10)
                lambda_group[s] ~ dgamma(0.1, 0.1)
            }
            # Individual level model: observations and process
            for(t in 1:T){
                for(s in 1:S){
                    for(j in 1:J){
                        for(i in 1:M){
                            z[i,j,s,t] ~ dbern(psi[j,s,t]) # Data augmentation variables
                            d[i,j,s,t] ~ dunif(0, B) # Distance is uniformly distributed
                            group_size[i,j,s,t] ~ dpois(lambda_group[s])  # Group size is Poisson
                            log(sigma[i,j,s,t]) <- alpha0[observer[i,j,s,t]] +  alpha1 * group_size[i,j,s,t]
                            mu[i,j,s,t] <- z[i,j,s,t]*exp(-d[i,j,s,t]*d[i,j,s,t]/(2*sigma[i,j,s,t]*sigma[i,j,s,t])) #p dep on dist class
                            # here using the half normal detection function (Buckland et al. 2001)
                            y[i,j,s,t] ~ dbern(mu[i,j,s,t])
                            zg[i,j,s,t]<- z[i,j,s,t]*(group_size[i,j,s,t] + 1)      # Number of individuals in that group
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
                lambda_grsize[s] <- lambda_group[s] + 1
                for(t in 1:T){
                    lambda_pop[s,t] <- exp(log_lambda_pop[s,t])
                    EN_region[s,t] <- lambda_grsize[s] * lambda_pop[s,t] * region_area[s]
                    #            for(j in 1:J){
                    #                realNg_transect[j,s,t] <- sum(z[,j,s,t]) # Realized transect-level estimates of number groups.
                    #                realNi_transect[j,s,t] <- sum(zg[,j,s,t])  # Total population size (all groups combined)
                    # NOTE: should not be assumed to be independent (i.e., do NOT sum across transects for population-level abundance estimate)
                    #            }
                    #            realNg_surveyregion[s,t] <- sum(realNg_transect[,s,t])
                    #            realDg_surveyregion[s,t] <- realNg_surveyregion[s,t]/region_area[s]
                    #            realNi_surveyregion[s,t] <- sum(realNi_transect[,s,t])
                    #            realDi_surveyregion[s,t] <- realNi_surveyregion[s,t]/region_area[s]
                }
            }
        }
        ",
        fill = TRUE, file = model_filepath)
}
