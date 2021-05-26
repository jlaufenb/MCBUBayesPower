

#' Write JAGS Model 2 to File
#'
#' @param model_filepath File path and name for JAGS model to be written
#'
#' @return Text file containing JAGS model is saved to \code{model_filepath}
#' @export
#'
write_model_2 <- function(model_filepath = "./models/model_2.txt"){
    cat("
        model{
            # single year, island-specific density and group size, random observer-effect on sigma model
            # Prior distributions for model parameters
            mu_obs ~ dnorm(0, 0.001)
            tau_obs <- pow(sd_obs,-2)
            sd_obs ~ dunif(0,10)
            for(k in 1:(K-1)){
                alpha0[k] ~ dnorm(mu_obs,tau_obs) # scale parameter, observer-specific intercepts
            }
            alpha0[K] <- 0
            alpha1 ~ dunif(-5,5) # scale parameter, group-size effect
            for(s in 1:S){
                log_lambda_pop[s] ~ dunif(0,10)
                lambda_group[s] ~ dgamma(0.1, 0.1)
            }
            # Individual level model: observations and process
            for(s in 1:S){
                for(j in 1:J){
                    for(i in 1:M){
                        z[i,j,s] ~ dbern(psi[j,s]) # Data augmentation variables
                        d[i,j,s] ~ dunif(0, B) # Distance is uniformly distributed
                        group_size[i,j,s] ~ dpois(lambda_group[s])  # Group size is Poisson
                        log(sigma[i,j,s]) <- alpha0[observer[i,j,s]] +  alpha1 * group_size[i,j,s]
                        mu[i,j,s] <- z[i,j,s]*exp(-d[i,j,s]*d[i,j,s]/(2*sigma[i,j,s]*sigma[i,j,s])) #p dep on dist class
                        # here using the half normal detection function (Buckland et al. 2001)
                        y[i,j,s] ~ dbern(mu[i,j,s])
                        zg[i,j,s]<- z[i,j,s]*(group_size[i,j,s] + 1)      # Number of individuals in that group
                    }
                }
            }
            for(s in 1:S){
                # Model for population size of groups
                for(j in 1:J){
                    Nj[j,s] ~ dpois(lambda_offset[j,s])
                    log(lambda_offset[j,s]) <- log_lambda_pop[s] + log(transect_area[j,s])
                    # psi is a derived parameter
                    psi[j,s] <- lambda_offset[j,s]/M
                }
            }
            # Derived quantities
            for(s in 1:S){
                lambda_pop[s] <- exp(log_lambda_pop[s])
                lambda_grsize[s] <- lambda_group[s] + 1
                EN_region[s] <- lambda_grsize[s] * lambda_pop[s] * region_area[s]
                for(j in 1:J){
                    realNg_transect[j,s] <- sum(z[,j,s]) # Realized transect-level estimates of number groups.
                    realNi_transect[j,s] <- sum(zg[,j,s])  # Total population size (all groups combined)
                    # NOTE: should not be assumed to be independent (i.e., do NOT sum across transects for population-level abundance estimate)
                }
                realNg_surveyregion[s] <- sum(realNg_transect[,s])
                realDg_surveyregion[s] <- realNg_surveyregion[s]/region_area[s]
                realNi_surveyregion[s] <- sum(realNi_transect[,s])
                realDi_surveyregion[s] <- realNi_surveyregion[s]/region_area[s]
            }
        }
        ",
        fill = TRUE, file = model_filepath)
}
