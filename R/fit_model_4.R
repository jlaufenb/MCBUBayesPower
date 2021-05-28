

#' Fit JAGS Model 4
#'
#' @param mcbu_data_file File path and name for MCBU distance data.
#' @param ni Total number of MCMC samples to simulate including burn-in samples. (Does not include adaptive samples)
#' @param nb Number of MCMC samples to discard as burn-in samples.
#' @param nt Integer specifying the thinning rate for MCMC samples.
#' @param nc Integer specifying the number of individual MCMC chains to simulate.
#' @param M Integer vector specifying the augmented number of MCBU groups per square kilometer.
#' @param save_output Logical value specifying whether to save JAGS output.
#' @param jags_output_file  File path and name used to save JAGS output.
#' @param use_parallel Logical value specifying whether to use parallel computing when simulating >1 MCMC chain. See help for \code{jags} function in \code{jagsUI} package for more details.
#'
#' @return  JAGS MCMC output from fitting an exponential growth model with individual birds as observation unit and observer random effect on sigma parameter to MCBU distance-sampling data from 2003 and 2018
#' @export
#'
fit_model4 <- function(mcbu_data_file = "./data/distdata_mcbu.csv", ni = 200, nb = 100, nt = 1, nc = 1, M = c(200,100),
                      save_output = FALSE, jags_output_file = NULL, use_parallel = TRUE){
    if(!"models" %in% list.files("./"))
        dir.create("./models")
    write_model_4()
    if(is.null(jags_output_file)){
        if(!"output" %in% list.files("./"))
            dir.create("./output")
        jags_output_file = "./output/out4.RData"
    }
    mcbu_data = read.csv(mcbu_data_file) # MCBU detections only, both survey years, data cleaned and ready for MCBU analysis
    exp_names = c("Region.Label","Area","Sample.Label","Effort","species","distance","size","obs","date","julian","year","bird_lat","bird_long","notes")
    if(!all(exp_names %in% colnames(mcbu_data))){
        stop(paste0("Column names missing from MCBU distance data. Expected column names are:\n\n",
                    paste0(exp_names, collapse = ", ")))
    }
    colnames(mcbu_data) = c("region_id","region_area","survey_number","transect_length","species","distance","group_size","observer",
                                 "date","julian","year","latitude","longitude","notes")
    indz = unlist(mapply(function(i,n){rep(i,each = n)}, i = seq_len(nrow(mcbu_data)),
                      n = as.integer(mcbu_data$group_size), SIMPLIFY = TRUE))
    mcbu_data = mcbu_data[indz,]
    mcbu_data$transect_id = paste0(mcbu_data$region_id,"_",mcbu_data$survey_number)
    mcbu_data$survey_number = as.character(mcbu_data$survey_number)
    mcbu_data$distancekm = mcbu_data$distance / 1000
    observer_levels = data.frame(observer = c("MAL","REG","SMM","JAJ","RMR","ARD","MND","BWR","MDR","DRR","SLW"),
                                 level = 11:1)
    mcbu_data$observer = as.factor(match(mcbu_data$observer, observer_levels$observer))
    data2003 = subset(mcbu_data, mcbu_data$year==2003)
    data2018 = subset(mcbu_data, mcbu_data$year==2018)
    split_data2003 = lapply(split(data2003, data2003$region_id, drop = TRUE),function(x)split(x, x$transect_id))
    split_data2018 = lapply(split(data2018, data2018$region_id, drop = TRUE),function(x)split(x, x$transect_id))
    datalist = list("y2003" = split_data2003, "y2018" = split_data2018)
    #
    yearstep = matrix(c(15,15),2,1)
    region_area = range(mcbu_data$region_area)
    K = nrow(observer_levels) + 1
    B = 0.1 # half-width of survey transects in km
    J = sapply(datalist, function(x)sapply(x,length))
    S = nrow(J)
    T = ncol(J)
    y = z = d = observer = array(0,c(max(M),max(J),S,T))
    z = d = observer = array(NA,c(max(M),max(J),S,T))
    transect_area = array(0,c(max(J),S,T))
    for(t in 1:T){
        for(s in 1:S){
            for(j in 1:J[s,t]){
                nobs = nrow(datalist[[t]][[s]][[j]])
                y[1:nobs,j,s,t] = z[1:nobs,j,s,t] = 1
                d[1:nobs,j,s,t] = datalist[[t]][[s]][[j]]$distancekm
                transect_area[j,s,t] = datalist[[t]][[s]][[j]]$transect_length[1] * B * 2
                observer[,j,s,t] = as.integer(datalist[[t]][[s]][[j]]$observer[1])
            }
        }
        z[,(min(J[,t])+1):max(J[,t]),1,t] = 0
        z[(min(M)+1):max(M),,1,t] = 0
        observer[,(min(J[,t])+1):max(J[,t]),1,t] = K
    }
    # format jags input data
    jags.data = list(M = max(M),
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
    params4 = c("mu_obs", "sd_obs", "alpha0", "r", "lambda_pop", "EN_region")
    inits = function(){list(mu_obs = -3, sd_obs = 1, r = 0, log_lambda0 = rep(4,S))}
    out4 = jagsUI::jags(jags.data, inits, params4, "./models/model_4.txt", n.thin = nt,
                 n.chains = nc, n.burnin = nb, n.iter = ni, parallel = use_parallel)
    if(save_output) save(out4, file = jags_output_file)
    return(out4)
}

