library(DoseFinding)
library(BB)
library(brglm)
library(numDeriv)
library(maxLik)
library(formula.tools)


###########Compute Hessian component of emax model####
Hess_e02 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1* sum(weight*pr*(1-pr))

  h
}

Hess_emax2 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1*sum(weight* pr*(1-pr)*(data$dose/(exp(led50)+data$dose))^2)

  h
}

Hess_led502 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1 * sum(weight * (pr*(1-pr)*(emax*data$dose*exp(led50)/(exp(led50)+data$dose)^2)^2+
                            (data$y-pr)*emax*data$dose*exp(led50)*(data$dose-exp(led50))/(exp(led50)+data$dose)^3))

  h
}

Hess_e0emax <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- -1* sum(weight* pr*(1-pr)*data$dose/(exp(led50)+data$dose))

  h
}

Hess_e0led50 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <- sum(weight* pr*(1-pr)*data$dose*emax/(exp(led50)+data$dose)^2*exp(led50))

  h
}

Hess_emaxled50 <- function(data,theta,weight){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  h <-  sum(weight* (pr*(1-pr)*data$dose^2*emax/(exp(led50)+data$dose)^3-(data$y-pr)*data$dose/(exp(led50)+data$dose)^2)*exp(led50))

  h
}

Comp_Hess <- function(data,theta,weight){
  theta <- as.numeric(theta)
  Hess <- matrix(0,3,3)
  Hess[1,1] <- Hess_e02(data,theta,weight)
  Hess[2,2] <- Hess_emax2(data,theta,weight)
  Hess[3,3] <- Hess_led502(data,theta,weight)
  Hess[1,2] <- Hess_e0emax(data,theta,weight)
  Hess[1,3] <- Hess_e0led50(data,theta,weight)
  Hess[2,3] <- Hess_emaxled50(data,theta,weight)
  Hess[3,1] <- Hess[1,3]
  Hess[3,2] <- Hess[2,3]
  Hess[2,1] <- Hess[1,2]

  Hess
}

Comp_Hess_deriv <- function(Hess,data,theta,weight){
  theta <- as.numeric(theta)
  Hess_deriv_e0 <- matrix(0,3,3)
  Hess_deriv_emax <- matrix(0,3,3)
  Hess_deriv_led50 <- matrix(0,3,3)

  grad_e02<- grad(Hess_e02,data=data,weight=weight,x=theta)
  grad_emax2<- grad(Hess_emax2,data=data,weight=weight,x=theta)
  grad_led502<- grad(Hess_led502,data=data,weight=weight,x=theta)
  grad_e0emax<- grad(Hess_e0emax,data=data,weight=weight,x=theta)
  grad_e0led50<- grad(Hess_e0led50,data=data,weight=weight,x=theta)
  grad_emaxled50<- grad(Hess_emaxled50,data=data,weight=weight,x=theta)

  Hess_deriv_e0[1,1] <- grad_e02[1]
  Hess_deriv_e0[2,2] <- grad_emax2[1]
  Hess_deriv_e0[3,3] <- grad_led502[1]
  Hess_deriv_e0[1,2] <- grad_e0emax[1]
  Hess_deriv_e0[1,3] <- grad_e0led50[1]
  Hess_deriv_e0[2,3] <- grad_emaxled50[1]
  Hess_deriv_e0[3,1] <- Hess_deriv_e0[1,3]
  Hess_deriv_e0[3,2] <- Hess_deriv_e0[2,3]
  Hess_deriv_e0[2,1] <- Hess_deriv_e0[1,2]

  Hess_deriv_emax[1,1] <- grad_e02[2]
  Hess_deriv_emax[2,2] <- grad_emax2[2]
  Hess_deriv_emax[3,3] <- grad_led502[2]
  Hess_deriv_emax[1,2] <- grad_e0emax[2]
  Hess_deriv_emax[1,3] <- grad_e0led50[2]
  Hess_deriv_emax[2,3] <- grad_emaxled50[2]
  Hess_deriv_emax[3,1] <- Hess_deriv_emax[1,3]
  Hess_deriv_emax[3,2] <- Hess_deriv_emax[2,3]
  Hess_deriv_emax[2,1] <- Hess_deriv_emax[1,2]

  Hess_deriv_led50[1,1] <- grad_e02[3]
  Hess_deriv_led50[2,2] <- grad_emax2[3]
  Hess_deriv_led50[3,3] <- grad_led502[3]
  Hess_deriv_led50[1,2] <- grad_e0emax[3]
  Hess_deriv_led50[1,3] <- grad_e0led50[3]
  Hess_deriv_led50[2,3] <- grad_emaxled50[3]
  Hess_deriv_led50[3,1] <- Hess_deriv_led50[1,3]
  Hess_deriv_led50[3,2] <- Hess_deriv_led50[2,3]
  Hess_deriv_led50[2,1] <- Hess_deriv_led50[1,2]

  return(list(Hess_grad_e0=Hess_deriv_e0,Hess_grad_emax=Hess_deriv_emax,Hess_grad_led50=Hess_deriv_led50))
}

Comp_I <- function(data,weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)




  I <- matrix(0,3,3)
  for(i in 1:length(data$y)){
    I_i <- matrix(0,3,3)
    I_i[1,] <- c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)
    I_i[2,] <- I_i[1,]*data$dose[i]/(exp(led50)+data$dose[i])
    I_i[3,] <- I_i[2,]*exp(led50)*(-emax)/(exp(led50)+data$dose[i])

    I_i <- I_i*weight[i]*pr[i]*(1-pr[i])

    I <- I+I_i
  }

  I
}
####compute theta with Firth likelihood


#' Jeffreys-penalized likelihood estimator via Newton–Raphson
#'
#' Maximizes the Jeffrey's prior-penalized log-likelihood for the binary Emax model
#' using \code{maxLik::maxNR}.
#'
#' @param data A data.frame (or list) with \code{y} (0/1) and \code{dose}.
#' @param weight Numeric vector of case weights.
#' @param theta Numeric(3) initial guess \code{c(e0, emax, led50)}.
#'
#' @details
#' Requires helpers \code{Comp_Hess()}, \code{Comp_I()}, and \code{Comp_Hess_deriv()}
#' that compute the (penalized) Hessian, information, and derivatives of the Hessian.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Estimated parameter vector.}
#'   \item{hessian}{Final Hessian returned by optimizer.}
#'   \item{vc}{Fisher-based variance/covariance.}
#' }
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#'  colnames(theta_true)<- c('e_0','emax','led_50')
#'  theta_true <- as.data.frame(theta_true)
#'  dose_set <- c(0,7.5,22.5,75,225)
#'  n=355
#'  data <-sim_data(theta_true,n,dose_set)
#'  res <- comp_theta_jeffrey(data=data )
#'  }
#' }
#' @seealso
#'  \code{\link[clinDR]{fitEmax}}
#' @export
comp_theta_jeffrey <- function(data=NULL,weight=NULL,theta=NULL){
  if(is.null(weight)){ weight <-rep(1,nrow(data))}
  theta_0 <-  clinDR::startEmax(y=data$y,dose=data$dose,binary = TRUE)
  theta_0 <- theta_0[c(3,2,1)]

  if(is.null(theta)){theta<-theta_0}

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  qr <- function(data,weight,theta){
    theta <- as.numeric(theta)
    #inout Hess

    Hess <- Comp_Hess(data,theta,weight)
    I <- Comp_I(data,weight,theta)
    logdetI <- ifelse(is.na(log(det(-Hess))),log(det(I)),log(det(-Hess)))

    # I <- Comp_I(data,theta,weight)
    # logdetI <- ifelse(is.na(log(det(I))),0,log(det(I)))
    e0 <- theta[1]
    emax <- theta[2]
    led50 <- theta[3]
    dose <- data$dose
    y <- data$y

    b <- -emax  #emax
    c <- 1/exp(led50) #ed50 related
    lambda <- 1
    a <- e0 - b #e0/emax
    evec <- a + b/(1 + (c * dose)^lambda)

    ll1 <- 0
    ll0 <- 0
    if (length(y == 1) > 0) {
      esub <- evec[y == 1]
      countsub <- weight[y == 1]
      ll1 <- sum(countsub * (plogis(esub, log.p = TRUE)))
    }
    if (length(y == 0) > 0) {
      esub <- evec[y == 0]
      countsub <- weight[y == 0]
      ll0 <- sum(countsub * (plogis(-esub, log.p = TRUE)))
    }
    return((ll1 + ll0+logdetI/2))
  }

  grr <- function(data,weight,theta){
    theta <- as.numeric(theta)

    e0 <- theta[1]
    emax <- theta[2]
    led50 <- theta[3]

    pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)


    #compute derivative of log detmerminant of hessian


    Hess <- Comp_Hess(data,theta,weight)
    grad_Hess <- Comp_Hess_deriv(Hess,data,theta,weight)

    I_inv <- MASS::ginv(-Hess)
    I_part_e0 <- (-1)*sum(diag(I_inv%*%grad_Hess$Hess_grad_e0))/2
    I_part_emax <- (-1)*sum(diag(I_inv%*%grad_Hess$Hess_grad_emax))/2
    I_part_led50 <- (-1)*sum(diag(I_inv%*%grad_Hess$Hess_grad_led50))/2


    q1 <- sum(weight*((data$y-pr)))+sum(I_part_e0)
    q2 <- sum(weight*((data$y-pr)*data$dose/(exp(led50)+data$dose)))+sum(I_part_emax)
    q3 <- sum(weight*((data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2))+sum(I_part_led50)

    c(q1,q2,q3)
  }

  # fit <- nlm(f = qr, p = as.numeric(c(e0,emax,ed50)), hessian = TRUE,
  #            data=data,weight=weight)

  #res<-spg(par=as.numeric(c(e0,emax,ed50)),fn=qr, lower=c(-Inf,-Inf,10^-1), control = list(maxit = 1e5),data=data,weight=weight)
  # res <- constrOptim(theta = as.numeric(c(e0,emax,led50)), f=qr,
  #                    grad= NULL,
  #                    ui=rbind(c(-1,1,0),  # the y-x > 0
  #                             c(0,0,1),   # the z > 0.001*maxdose
  #                             c(0,0,-1)), # the z<1.2*maxdose
  #                    ci=c(0,log(0.001*max(data$dose)),-log(1.2*max(data$dose))),data=data,weight=weight)
  #
  res <- maxLik::maxNR(fn=qr, grad=NULL,start =as.numeric(c(e0,emax,led50)),data=data,weight=weight )

  #H<-hessian(qr,x= as.numeric(res$par),data=data,weight=weight)
  # return(res)

  I <- Comp_I(data=data,weight=weight,theta=res$estimate)
  vc <- try(MASS::ginv(-res$hessian))
  if (inherits(vc,'try-error')){vc=matrix(NA,3,3)}
  return(list(par=res$estimate,hessian=res$hessian,vc=vc))
}






#' @title Simulate dataset for testing three bias-corrected meethods
#' @description FUNCTION_DESCRIPTION
#' @param theta True parameters of Emax model. The order of the variables is (E0,log(ED50),Emax).
#' @param n number of observations.
#' @param dose_set A vector indicate the dose set for the dose-response relationship.
#' @return A dataframe of simulated dataset.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#' colnames(theta_true)<- c('e_0','emax','led_50')
#' theta_true <- as.data.frame(theta_true)
#' dose_set <- c(0,7.5,22.5,75,225)
#' n=355

#' data <-sim_data(theta_true,n,dose_set)

#'  }
#' }
#' @rdname sim_data
#' @export

sim_data <- function(theta,n,dose_set){
  dose <- rep(dose_set,each=round(n/length(dose_set)))
  pi <- plogis(theta$e_0+theta$emax*dose/(exp(theta$led_50)+dose))
  y <- rbinom(n,1,pi)
  id <- seq(1,length(y),1)
  data <- cbind(id,y,dose)
  return(as.data.frame(data))
}
