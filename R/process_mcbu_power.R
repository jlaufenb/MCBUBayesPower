

#' Process and Summarize MCBU Power Simulation Output
#'
#' @param true_pctdecl True percent decline over 10 years used in simulations.
#' @param jags_filepath File path to folder containing saved JAGS output from power simulations.
#' @param save_figures Logical value specifying whether to save summary figures.
#' @param figure_filepath File path to folder to which figures are saved. Default creates a "figures" folder along side the folder containing saved JAGS output from power simulations.
#' @param save_summary Logical value specifying whether to save summary R objects.
#' @param summary_filepath File path to folder to which summary R objects are saved. Default creates a "summary" folder along side the folder containing saved JAGS output from power simulations.
#'
#' @return None
#' @export
#'
process_mcbu_power <- function(true_pctdecl = 25, output_filepath = NULL,
                               save_figures = TRUE, figure_filepath = NULL,
                               save_summary = TRUE, summary_filepath = NULL){

    if(is.null(output_filepath))stop("Must provide file path to output_filepath argument for the folder containing RData files with JAGS output.")
    filenames = list.files(output_filepath, full.names = TRUE)
    nreps = length(filenames)
    true_lambda = (1 - (true_pctdecl/100))^(1/10)
    true_r = log(true_lambda)
    cv = bias = cover_95 = onetail_95 = cover_90 = onetail_90 = rep(0,nreps)
    pctdecl_ests = matrix(0,nreps,length(summary.fn()))
    colnames(pctdecl_ests) = summary.fn()
    print(system.time({
        for(rep in 1:nreps){
            cat("doing", rep, " of ", nreps, " reps\n\n")
            fname = filenames[rep]
            jags_rep = get(load(fname))
            rm(list = paste0("mcmc",substr(fname, regexpr("mcmc_",fname)[1] + 5, regexpr(".RData",fname)[1]-1)))
            r_post = jags_rep$sims.list$r
            r_est = summary.fn(r_post)
            lambda_post = exp(r_post)
            lambda_est = summary.fn(lambda_post)
            pctdecl_post = (1 - lambda_post^10) * 100
            pctdecl_est = summary.fn(pctdecl_post)
            pctdecl_ests[rep,] = pctdecl_est
            nsamples = length(r_post)
            onetail_95[rep] = sum(pctdecl_post < true_pctdecl)/nsamples < 0.05 # one-tailed test for at least 25% decline at 0.95 prob
            onetail_90[rep] = sum(pctdecl_post < true_pctdecl)/nsamples < 0.1 # one-tailed test for at least 25% decline at 0.90 prob
            cv[rep] = as.numeric(abs(pctdecl_est["sd"]/pctdecl_est["mean"]))
            bias[rep] = as.numeric(pctdecl_est["mean"] - true_pctdecl)
            lims = pctdecl_est[c("2.5%","97.5%")]
            cover_95[rep] = lims[1] < true_pctdecl & lims[2] > true_pctdecl
            lims = pctdecl_est[c("5%","95%")]
            cover_90[rep] = lims[1] < true_pctdecl & lims[2] > true_pctdecl
        }
    }))
    #
    lty_indz <- -(cover_95 - 1) + 1
    sq_err = bias^2
    sq_err_sum = summary.fn(sq_err)
    cv_sum = summary.fn(cv)
    bias_sum = summary.fn(bias)
    rmse = sqrt(mean(sq_err))
    (ci_cover_95 = sum(cover_95)/nreps)
    (onetail_power_95 = 1 - (sum(onetail_95)/nreps))
    (ci_cover_90 = sum(cover_90)/nreps)
    (onetail_power_90 = 1 - (sum(onetail_90)/nreps))
    #
    if(save_summary){
        if(is.null(summary_filepath)){
            indz = gregexpr("/",output_filepath)[[1]]
            nindz = length(indz)
            fig_dir = paste0(substr(output_filepath,1,indz[nindz-1]),"summary")
            if(!dir.exists(fig_dir))dir.create(fig_dir)
        }
        save(sq_err, cv, bias, rmse, cover_95, cover_90, ci_cover_95, ci_cover_90, onetail_power_95, onetail_power_90,
                         r_post, lambda_post, pctdecl_post,pctdecl_est,
                         file = paste0(summary_filepath, "/summary.RData"))
    }
    if(save_figures){
        if(is.null(figure_filepath)){
            indz = gregexpr("/",output_filepath)[[1]]
            nindz = length(indz)
            fig_dir = paste0(substr(output_filepath,1,indz[nindz-1]),"figures")
            if(!dir.exists(fig_dir))dir.create(fig_dir)
        }

        png(paste0(figure_filepath, "/ci_cover.png"), height = 9, width = 6.5, units = "in", res = 192)
        plot(1, type = "n", ylim = c(0, nreps), xlim = c(0,50), ylab = "Iteration", xlab = "Percent decline", main = "", frame.plot = FALSE)
        segments(pctdecl_ests[,"2.5%"], 1:nreps, pctdecl_ests[,"97.5%"], 1:nreps,
                 col = c("blue","red")[lty_indz], lty = c(1,3)[lty_indz], lwd = c(1,1.5)[lty_indz])
        abline(v = true_pctdecl)
        abline(v = mean(pctdecl_ests[,"mean"]), lty = 2)
        points(pctdecl_ests[,"mean"], 1:nreps, pch = 17, cex = 0.5)
        dev.off()
        #
        png(paste0(figure_filepath, "/cv_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
        hist(cv, xlim = c(0,0.4), breaks = seq(0,0.4,0.005), xlab = "Coefficient of variation", main = "")
        abline(v = cv_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
        dev.off()
        #
        png(paste0(figure_filepath, "/bias_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
        hist(bias, xlim = c(-15,15), breaks = seq(-15, 15, 1), xlab = "Bias", main = "")
        abline(v = bias_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
        dev.off()
        #
        png(paste0(figure_filepath, "/sq_err_hist.png"), height = 9, width = 6.5, units = "in", res = 192)
        hist(sq_err, xlim = c(0,150), breaks = seq(0, 150, 5), xlab = "Squared error", main = "")
        abline(v = sq_err_sum[c("mean","2.5%","97.5%")], lty = c(1,2,2), lwd = 2, col = "blue")
        dev.off()
    }
}
