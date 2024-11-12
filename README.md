# tvvarnets
`tvvarnets.R` implements time-varying group lasso to estimate a uniform time-varying *Granger causality network* and time-varying CLIME to estimate a uniform time-varying *partial correlation network* based on a time-varying VAR model. 

See

> Chen, J., Li, D., Li, Y. and Linton, O., 2025. Estimating time-varying networks for high-dimensional time series.

## main 
*main_net_fac.R*:		time-varying network 

*main_LASSO_fac.R*:	constant network with factors

*main_LASSO.R*: 		constant network

* create a results folder to save results
* check working directory, and current version is for running on a cluster.
* dyn.load("code/rawfit_wglasso.so", local=TRUE) for Linux
* dyn.load("code/rawfit_wglasso.dll", local=TRUE) for Windows

## DGP

*simTOP.R*: Toeplitz matrix

*simUP.R*: upper triangular matrix

*simDiag.R*: diagonal matrix

*simUPFac.R* diagonal matrix

## C code for weighted group LASSO with extremely large design matrix 

rawfit_wglasso.c
