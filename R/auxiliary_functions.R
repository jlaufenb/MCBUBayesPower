#' Create Folders
#'
#' @param dir Location to create new 'data', 'models', and 'output' folders.
#' @param moreFolders Character vector of additional folder names to be created.
#'
#'@export

create_folders <- function(dir = NULL, moreFolders = NULL) {

    ## folder location and names
    if(is.null(dir))root.dir = getwd()
    folder_names = c("data","models","output",moreFolders)

    ## check for existing folders
    folder_names = folder_names[!(folder_names %in% list.files(root.dir))]

    ## create new folders
    f = paste0(root.dir, "/", folder_names)
    lapply(f, dir.create)
}



#' Calculate Mode
#'
#' @param x vector of numeric values
#'
#' @return numeric value of mode as estimated by the base \code{density} function in R with \code{na.rm = TRUE}
#' @export
#'
Mode <- function(x) {
    d = density(x,na.rm=TRUE)
    d$x[which.max(d$y)]
}

#' Calculate Simple Summary Statistics
#'
#' @param x numeric vector or matrix
#'
#' @return numeric vector or matrix when x is supplied. Returns labels when x is NULL.
#' @export
#'
summary.fn	<- function(x=NULL){
    cn = c("mean","sd","mode","0%","2.5%","5%","25%","50%","75%","95%","97.5%","100%")
    if(!is.null(x)){
        if(any(is.na(x))){
            cat("\n\nWarning: x contains missing values\n\n")
        } #if
        if(is.matrix(x)){
            y = array(0L,c(dim(x)[2],10))
            for(i in 1:dim(x)[2]){
                y[i,] = c(mean(x[,i],na.rm=TRUE), sd(x[,i],na.rm=TRUE), Mode(x[,i]),
                          quantile(x[,i], probs = c(0,0.025,0.05,0.25,0.50,0.75,0.95,0.975,1), na.rm=TRUE))
                dimnames(y)[[2]] = cn
                dimnames(y)[[1]] = dimnames(x)[[2]]
            } #i
        }else{
            y = c(mean(x,na.rm=TRUE), sd(x,na.rm=TRUE), Mode(x),
                  quantile(x, probs = c(0,0.025,0.05,0.25,0.50,0.75,0.95,0.975,1), na.rm=TRUE))
            names(y) = cn
        } #ifelse
        return(y)
    }
    if(is.null(x)){return(cn)}
}






#' Extract MCMC Samples for Expected N
#'
#' @param mod_output MCMC output object of class \code{jagsUI} from either model.
#'
#' @return Array containing MCMC samples for island-specific and total abundance for surveyed years (2003 and 2018) and projected abundances for 2028.
#'
getEN_mcmc <- function(mod_output = NULL){
    EN_mcmc <- array(0, dim = c(mod_output$mcmc.info$n.samples,3,3), dimnames = list(NULL, c("hall","stma","total"),c(2003, 2018, 2028)))
    EN_mcmc[,1:2,1:2] <- mod_output$sims.list$EN_region
    EN_mcmc[,3,1:2] <- t(apply(EN_mcmc[,1:2,1:2],1,colSums))
    EN_mcmc[,1:3,3] <- EN_mcmc[,1:3,2] * exp(out4$sims.list$r)^10
    return(EN_mcmc)
}
