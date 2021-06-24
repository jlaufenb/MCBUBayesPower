




#' Calculate Distance Detection Probabilities
#'
#' @param dists Numeric vector of distances from transect line.
#' @param sigma Spatial-scale parameter for Gaussian (half-normal) distance detection function.
#'
#' @return Numeric vector of detection probabilities.
#' @export
#'


distdetfn <- function(dists = NULL,sigma = NULL){exp(-(dists * dists) / (2 * sigma^2))}




###########################################################################################################################


#' Plot Distance Detection Function
#'
#' Function plots distance detection functions for individual models (i.e., input provided for 1 model argument)
#' or functions from both models (i.e., input provided for both model argument) overlaid on the sample plot.
#'
#' @param mod3_output MCMC output object of class \code{jagsUI} for model with groups as observation unit.
#' @param mod4_output MCMC output object of class \code{jagsUI} for model with individuals as observation unit.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


detfn_plot <-  function(mod3_output = NULL, mod4_output = NULL, save_plot = FALSE, file_path = NULL, w = 9, h = 6.5){
    output_list = list(mod3 = mod3_output, mod4 = mod4_output)
    null_test = sapply(output_list, is.null)
    if(all(null_test)) stop("Warning: Must provide output from at least one model")
    output_list = output_list[!null_test]
    if(!any(sapply(output_list, class) %in% "jagsUI"))
        stop("Warning: Input data objects must be of class 'jagsUI'.")
    sig_list = vector("list",length(output_list))
    for(i in 1:length(sig_list)){
        sig_list[[i]] = c(exp(output_list[[i]]$mean$mu_obs) * 1000,
                          exp(output_list[[i]]$mean$mu_obs + 1.96 * output_list[[i]]$mean$sd_obs)*1000,
                          exp(output_list[[i]]$mean$mu_obs - 1.96 * output_list[[i]]$mean$sd_obs)*1000,
                          exp(output_list[[i]]$mean$alpha0[1:11]) * 1000)
    }
    dists = seq(0,250,0.5)
    if(length(output_list)==2){
        if(is.null(file_path)) file_path = paste0(getwd(),"/detfn_model_compare")
        detfn_plot = function(){
            plot(1,type="n",xlim = c(0,250), ylim = c(0,1), frame.plot = FALSE, xlab = "Distance (km)", ylab = "Detection probability")
            polygon(c(dists,rev(dists)), c(distdetfn(dists, sig_list[[1]][3]),distdetfn(rev(dists),sig_list[[1]][2])), col = rgb(0.53, 0.81, 0.92, 0.5), lwd = 1.5, lty = 2)
            lines(dists,distdetfn(dists,sig_list[[1]][1]), lwd = 1.5, lty = 2)
            polygon(c(dists,rev(dists)), c(distdetfn(dists, sig_list[[2]][3]),distdetfn(rev(dists),sig_list[[2]][2])), col = rgb(0.53, 0.81, 0.92, 0.5), lwd = 1.5, lty = 1)
            lines(dists,distdetfn(dists,sig_list[[2]][1]), lwd = 1.5, lty = 1)
        }
    }else{
        if(is.null(file_path)) file_path = paste0(getwd(),"/detfn_",names(output_list))
        diffs <- ((output_list[[1]]$mean$alpha0 - output_list[[1]]$mean$mu_obs)/output_list[[1]]$mean$mu_obs) * -100
        detfn_plot = function(){
            plot(1,type="n",xlim = c(0,250), ylim = c(0,1), frame.plot = FALSE, xlab = "Distance (km)", ylab = "Detection probability")
            polygon(c(dists,rev(dists)), c(distdetfn(dists, sig_list[[1]][3]),distdetfn(rev(dists),sig_list[[1]][2])), col = "skyblue", lwd = 1.5, lty = 1)
            for(i in 4:length(sig_list[[1]])){
                lines(dists,distdetfn(dists,sig_list[[1]][i]), lwd = 1, lty = 3)
            }
            lines(dists,distdetfn(dists,sig_list[[1]][1]), lwd = 2)
            text((2:12)*20, 1.02, labels = paste0(sort(round(diffs[1:11],0)), "%"), cex = 0.75, adj = 0.5)
        }
    }
    if(save_plot){
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        detfn_plot()
        dev.off()
    }else{
        detfn_plot()
    }
}



###########################################################################################################################


#' Plot Posterior of Density Difference Between Islands
#'
#' @param mod_output MCMC output object of class \code{jagsUI} from either model.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


density_diff_plot <- function(mod_output = NULL, save_plot = FALSE, file_path = NULL, w = 6.5, h = 6.5){
    if(!class(mod_output) %in% "jagsUI")stop("Warning: Input data object must be of class 'jagsUI'.")
    lambda_pops <- mod_output$sims.list$lambda_pop
    lambda_pops_diffs <- lambda_pops_pctdiffs <- matrix(0,dim(lambda_pops)[1],2)
    for(i in 1:dim(lambda_pops)[1]){
        lambda_pops_diffs[i,] <- (lambda_pops[i,1,] - lambda_pops[i,2,])
        lambda_pops_pctdiffs[i,] <- (lambda_pops_diffs[i,] / lambda_pops[i,2,]) * 100
    }
    pctdiffs_summ = summary.fn(lambda_pops_pctdiffs[,1])
    modname = ifelse("lambda_group" %in% names(mod_output$mean), "mod3","mod4")
    diff_plot = function(){
        plot(1, type = "n", main = "", xlim = c(pctdiffs_summ["0%"],pctdiffs_summ["100%"]), ylim = c(0,0.08), xlab = "Percent", frame.plot = FALSE, ylab = "Density")
        x <- seq(pctdiffs_summ["0%"],pctdiffs_summ["100%"],0.1)
        y <- dnorm(x,pctdiffs_summ["mean"],pctdiffs_summ["sd"])
        polygon(x, y, col = rgb(0.53, 0.81, 0.92), lwd = 1, lty = 1)
        lines(rep(pctdiffs_summ["mean"],2),c(0,max(y)))
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/density_diff_", modname)
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        diff_plot()
        dev.off()
    }else{
        diff_plot()
    }

}








###########################################################################################################################




#' Plot Bias Histogram
#'
#' @param bias Numeric vector containing values of relative bias between posterior means and true value of the percent population decline per 10 years across all simulated data sets.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


bias_plot <- function(bias = NULL, save_plot = FALSE, file_path = NULL, w = 6.5, h = 6.5){
    bias_hist = function(){
        hist(bias, xlim = c(-15,15), breaks = seq(-15, 15, 1), xlab = "Bias", main = "")
        abline(v = bias_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/bias_histogram")
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        bias_hist()
        dev.off()
    }else{
        bias_hist()
    }
}








###########################################################################################################################



#' Plot Coefficient of Variation Histogram
#'
#' @param cv Numeric vector containing coefficient of variation (CV) estimates calculated from posterior distributions for the percent population decline per 10 years across all simulated data sets.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


cv_plot <- function(cv = NULL, save_plot = FALSE, file_path = NULL, w = 6.5, h = 6.5){
    cv_hist = function(){
        hist(cv, xlim = c(0,0.4), breaks = seq(0,0.4,0.005), xlab = "Coefficient of variation", main = "")
        abline(v = cv_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/cv_histogram")
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        cv_hist()
        dev.off()
    }else{
        cv_hist()
    }
}





###########################################################################################################################



#' Plot Bayesian Credible Intervals for All Simulations
#'
#' @param pctdecl_ests Matrix containing summary statistics calculated from posterior distributions of the percent population decline per 10 years for all simulated data sets.
#' @param cover_95 Numeric vector containing posterior distributions means of the percent population decline per 10 years for all simulated data sets.
#' @param true_pctdecl Numeric value specifying true percent population decline per 10 years used for simulations.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


ci_cover_plot <- function(pctdecl_ests = NULL, cover_95 = NULL, true_pctdecl = 25, save_plot = FALSE, file_path = NULL, w = 6.5, h = 9){
    ci_cover = function(){
        lty_indz <- -(cover_95 - 1) + 1
        nreps = length(pctdecl_ests)
        plot(1, type = "n", ylim = c(0, nreps), xlim = c(0,50), ylab = "Iteration", xlab = "Percent decline", main = "", frame.plot = FALSE)
        segments(pctdecl_ests[,"2.5%"], 1:nreps, pctdecl_ests[,"97.5%"], 1:nreps,
                 col = c("blue","red")[lty_indz], lty = c(1,3)[lty_indz], lwd = c(1,1.5)[lty_indz])
        points(pctdecl_ests[,"mean"], 1:nreps, pch = 17, cex = 0.5)
        abline(v = true_pctdecl, col = "pink", lwd = 1.5)
        abline(v = mean(pctdecl_ests[,"mean"]), lty = 2, col = "pink", lwd = 1.5)
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/ci_cover_plot")
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        ci_cover()
        dev.off()
    }else{
        ci_cover()
    }
}




###########################################################################################################################





#' Plot Squared Error Histogram
#'
#' @param sq_err Numeric vector containing squared error estimates calculated from posterior distributions for the percent population decline per 10 years across all simulated data sets.
#' @param save_plot Logical value indicating whether to save plot.
#' @param file_path File path including user-provided file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


sq_err_plot <- function(sq_err = NULL, save_plot = FALSE, file_path = NULL, w = 6.5, h = 6.5){
    sq_err_hist = function(){
        hist(sq_err, xlim = c(0,150), breaks = seq(0, 150, 5), xlab = "Squared error", main = "")
        abline(v = sq_err_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/sq_err_histogram")
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        sq_err_hist()
        dev.off()
    }else{
        sq_err_hist()
    }
}








###########################################################################################################################


#' Plot Observer-Specific Detection Functions and Observations
#'
#' @param distdata Dataframe containing raw MCBU distance-sampling data.
#' @param mod_output MCMC output object of class \code{jagsUI} from either model.
#' @param save_plots Logical value indicating whether to save plots.
#' @param file_path File path including user-provided file name. Observer indexes will automatically be appended to file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'


obs_detfn_plot <- function(distdata = NULL, mod_output = NULL, save_plots = FALSE, file_path = NULL, w = 6.5, h = 6.5){
    dists <- seq(0,250,0.5)
    sigmas <- c(exp(mod_output$mean$mu_obs) * 1000,
                 exp(mod_output$mean$mu_obs + 1.96 * mod_output$mean$sd_obs)*1000,
                 exp(mod_output$mean$mu_obs - 1.96 * mod_output$mean$sd_obs)*1000,
                 exp(mod_output$mean$alpha0[1:11]) * 1000)

    observer_levels = data.frame(observer = c("MAL", "REG", "SMM", "JAJ", "RMR", "ARD", "MND", "BWR", "MDR", "DRR", "SLW"), level = 11:1)
    distdata$obs = as.factor(match(distdata$obs, observer_levels$observer))
    by_obs <- split(distdata, distdata$obs)
    obs_detfn = function(){
        par(mar=c(5, 4, 4, 6) + 0.1)
        plot(1, type = "n", xlim = c(0,200), ylim = c(0,nrow(tmp)), frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")
        points(tmp$distance, 1:nrow(tmp), cex = tmp$size/6, pch = 16)
        abline(v = 100, lty = 3)
        mtext("Observation", side = 2, line = 2.5)
        axis(2)
        par(new = TRUE)
        plot(dists[!dists>200],distdetfn(dists[!dists>200],sigmas[i+3]), ylim = c(0,1), type="l", axes = FALSE, xlab = "", ylab = "", lwd = 2)
        mtext("Detection probability", side = 4, line = 2.5)
        axis(4, ylim = c(0,1))
        mtext("Distance (m)", side = 1, line = 2.5)
        axis(1)
        lines(dists[!dists>200],distdetfn(dists[!dists>200],sigmas[1]), lwd = 2, lty = 3)
    }
    if(save_plots){
        if(is.null(file_path)) file_path = paste0(getwd(),"/observer_detfn")
        for(i in 1:length(by_obs)){
            tmp <- by_obs[[i]]
            png(paste0(file_path,"_",i,".png"), width = w, height = h, units = "in", res = 192)
            obs_detfn()
            dev.off()
        }
    }else{
        for(i in 1:length(by_obs)){
            tmp <- by_obs[[i]]
            obs_detfn()
            readline(prompt="Hit <Enter> to advance to next plot.\n")
        }
    }
}



###########################################################################################################################



#' Plot Posterior Sample Histograms for Expected Abundance
#'
#' @param mod_output MCMC output object of class \code{jagsUI} from either model.
#' @param estimate Character vector specifying one or more island estimates to plot. Accepted values are "total", "hall", and "stma" and default is "total".
#' @param year Numeric vector specifying years for which island estimates are plotted. Accepted values are  2003, 2018, and 2028 and the default is all 3.
#' @param save_plots Logical value indicating whether to save plots.
#' @param file_path File path including user-provided file name. Observer indexes will automatically be appended to file name.
#' @param w Numeric value specifying figure width in inches.
#' @param h Numeric value specifying figure height in inches.
#'
#' @return
#' @export
#'

EN_isl_plot <- function(mod_output = NULL, estimate = "total", year = NULL, save_plot = FALSE, file_path = NULL, w = 9, h = 6.5){
    ests = c("hall","stma","total")
    Ests = c("Hall Island","St Matthew Island","Total island")
    yrs = c(2003, 2018, 2028)
    if(any(!(estimate %in% ests)))stop("Warning: Accepted values for the estimate argument are 'total', 'hall', and 'stma'")
    if(is.null(year)){
        year = yrs
    }else{
        if(any(!(year %in% yrs)))stop("Warning: Only 2003, 2018, and 2028 are currently accepted values for the year argument")
    }
    if(length(estimate) > 1 & length(year) > 1)stop("Warning: Cannot plot multiple estimate types (e.g., estimate = c('hall','stma')) and mutliple years (e.g., year = c(2003, 2018))")
    rows = which(ests %in% estimate)
    cols = which(yrs %in% year)
    EN_mcmc = getEN_mcmc(mod_output)
    EN_stats = apply(EN_mcmc,c(2,3),summary.fn)
    indz = list(rows,cols)
    nindz = sapply(indz,length)
    nhists = max(unlist(nindz))
    EN_mcmc = EN_mcmc[,indz[[1]], indz[[2]]]
    EN_stats = as.data.frame(EN_stats[,indz[[1]], indz[[2]]])
    xlims = c(0,35000)
    brks = seq(0,35000,100)
    if(all(ests[rows] %in% "hall")){
        xlims = c(0,2500)
        brks = seq(0,2500,25)
    }
    if(!("hall" %in% ests[rows])){
        xlims = c(5000,35000)
        brks = seq(5000,35000,250)
    }
    EN_plot = function(){
        plot(1, type = "n", xlim = xlims, ylim = c(0,30000), frame.plot = FALSE, ylab = "Frequency", xlab = "Abundance")
        for(i in 1:nhists){
            hist(EN_mcmc[,i], breaks = brks, col = rgb(211/255,211/255,211/255,0.75), add = T)
        }
        segments(x0 = unlist(EN_stats["mean",]), y0 = rep(0,3), x1 = unlist(EN_stats["mean",]), y1 = rep(25000,3), lty = 1, lwd = 2, col = "red")
        segments(x0 = unlist(EN_stats[c("2.5%","97.5%"),]), y0 = rep(0,3), x1 = unlist(EN_stats[c("2.5%","97.5%"),]), y1 = rep(25000,3), lty = 2, lwd = 2, col = "blue")
        text(unlist(EN_stats["mean",]), 28000, paste0(yrs[cols], " ", Ests[rows], "\n N-hat = ", round(unlist(EN_stats["mean",]),0)), cex=0.4)
        # text(unlist(EN_stats["mean",]), 28000, paste0("N-hat = \n", round(unlist(EN_stats["mean",]),0)), cex=0.5)
    }
    if(save_plot){
        if(is.null(file_path)) file_path = paste0(getwd(),"/output/EN_posteriors_",
                                                  paste0(ests[rows], collapse = "-"), "_",
                                                  paste0(yrs[cols], collapse = "-"))
        png(paste0(file_path,".png"), width = w, height = h, units = "in", res = 192)
        EN_plot()
        dev.off()
    }else{
        EN_plot()
    }
}

