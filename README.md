#Jeffreys-Prior Penalized Binary Emax Model Package

In this document, we illustrate the main features of the `jpemax` R package through examples. Additional information on the statistical methodology and computational details are provided in the accompanying documentation and research articles.

## Cite the package

The package applies methods introduced in the [paper]:

Zhang, J., Pradhan, V. and Zhao, Y., 2025. Evaluating Bias Reduction Methods in Binary Emax Model for Reliable Dose-Response Estimation.

## Install

Open the R console and run the following command to install the package from source:

```r
install.packages("devtools") # When you have not installed devtools package
devtools::install_github("Celaeno1017/jpemax")
```

## Tutorial

First, load the R package.

```r
library(ememax)
```

To illustrate the main features of the R package `jpemax`, let's first generate some data. We have built in a few functions directly into the R package for this purpose.

```r
theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
 colnames(theta_true)<- c('e_0','emax','led_50')
 theta_true <- as.data.frame(theta_true)
 dose_set <- c(0,7.5,22.5,75,225)
 n=355
 data <-sim_data(theta_true,n,dose_set)
```

To fit the penalized Emax model with Jeffreys prior, we use the comp_theta_jeffrey function which implements the proposed methodology.

```r
res <- comp_theta_jeffrey(data=data )
```
Key parameters include:
- `data`: A data.frame (or list) with y (0/1) and dose.
- `weight`: Numeric vector of case weights.
- `theta`: Numeric(3) initial value of c(e0, emax, log(ed50)).
- 
The result will contain the following values:
- `par`:	the final fitted parameters of Emax model

- `hessian`:	the final fitted Hessian returned by optimizer.

- `vc`: Fisher information-based variance/covariance.
