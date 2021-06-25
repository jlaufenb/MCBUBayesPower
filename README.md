# USFWS Disclaimer
The United States Fish and Wildlife Service (FWS) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. FWS has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by FWS. The FWS seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by FWS or the United States Government.

# MCBUBayesPower
MCBUBayesPower is a custom package created to fit hierarchical population models to distance-sampling
data using Bayesian methods to obtain parameter estimates, to fit models to simulated data sets based
on prospective survey designs, and compile simulation results to evaluate statistical power among different
survey designs. 

# Instructions

To install and load the package:  

`if (!require("devtools")) install.packages("devtools")`  
`devtools::install_github("jlaufenb/MCBUBayesPower", build_vignettes = TRUE)`

`library(MCBUBayesPower)`

Further instructions describing the basic workflow for using `MCBUBayesPower` is provided as a vignette in the package help documentation:

help(package = "MCBUBayesPower")
