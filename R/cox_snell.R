comp_Q_e0 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(-1*weight^2* pr*(1-pr)*data$dose/(exp(led50)+data$dose)^2*exp(led50))
  Q[3,3] <- sum(weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax*data$dose*exp(led50)/(exp(led50)+data$dose)^3)
  Q[3,2] <- Q[2,3]

  Q
}

comp_Q_emax <- function(data,weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(-1*weight^2* pr*(1-pr)*data$dose^2/(exp(led50)+data$dose)^3*exp(led50))
  Q[3,3] <- sum(weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax*data$dose^2*exp(led50)/(exp(led50)+data$dose)^4)
  Q[3,2] <- Q[2,3]

  Q
}

comp_Q_led50 <- function(data,weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)

  Q <- matrix(0,3,3)
  Q[2,3] <- sum(weight^2* pr*(1-pr)*data$dose^2*emax/(exp(led50)+data$dose)^4*exp(led50)^2)
  Q[3,3] <- sum(-1*weight^2* pr*(1-pr)*(data$dose-exp(led50))*emax^2*data$dose^2*exp(led50)^2/(exp(led50)+data$dose)^5)
  Q[3,2] <- Q[2,3]

  Q
}

Comp_I_score <- function(data,weight, theta){
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


comp_K_e0 <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(2*pr[i]-1)*pr[i]*(1-pr[i])

    I_i[2,] <- I_i[2,]-c(0,0,-pr[i]*(1-pr[i])*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)
    I_i[3,] <-  I_i[3,]-c(0,-pr[i]*(1-pr[i])*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2,pr[i]*(1-pr[i])*data$dose[i]*emax*exp(led50)*(data$dose[i]-exp(led50))/(exp(led50)+data$dose[i])^3)
    I <- I+I_i
  }

  I

}

comp_K_emax <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(2*pr[i]-1)*pr[i]*(1-pr[i])*data$dose[i]/(exp(led50)+data$dose[i])

    I_i <- I_i - pr[i]*(1-pr[i])*( c(0,0,-data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)%*%t(c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2))+
                                     c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)%*%t(c(0,0,-data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)))


    I_i[2,] <- I_i[2,]-c(0,0,-pr[i]*(1-pr[i])*data$dose[i]^2*exp(led50)/(exp(led50)+data$dose[i])^3)
    I_i[3,] <-  I_i[3,]-c(0,-pr[i]*(1-pr[i])*data$dose[i]^2*exp(led50)/(exp(led50)+data$dose[i])^3,pr[i]*(1-pr[i])*data$dose[i]^2*emax*exp(led50)*(data$dose[i]-exp(led50))/(exp(led50)+data$dose[i])^4)

    I <- I+I_i
  }

  I

}


comp_K_led50 <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(2*pr[i]-1)*pr[i]*(1-pr[i])*(-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)

    I_i <- I_i - pr[i]*(1-pr[i])* (c(0,-data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2,-data$dose[i]*emax*exp(led50)*(data$dose[i]-exp(led50))/(exp(led50)+data$dose[i])^3)%*%
      t(c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2))+
        c(1,data$dose[i]/(exp(led50)+data$dose[i]),-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)%*%
        t(c(0,-data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2,-data$dose[i]*emax*exp(led50)*(data$dose[i]-exp(led50))/(exp(led50)+data$dose[i])^3)))

    I_i[2,] <- I_i[2,]-c(0,0,pr[i]*(1-pr[i])*data$dose[i]^2*exp(led50)^2*emax/(exp(led50)+data$dose[i])^4)
    I_i[3,] <-  I_i[3,]-c(0,pr[i]*(1-pr[i])*data$dose[i]^2*exp(led50)^2*emax/(exp(led50)+data$dose[i])^4,-pr[i]*(1-pr[i])*data$dose[i]^2*emax^2*exp(led50)^2*(data$dose[i]-exp(led50))/(exp(led50)+data$dose[i])^5)


    I <- I+I_i
  }

  I

}


comp_bias <- function(data, weight,theta){

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  Q_e0 <- comp_Q_e0(data, weight,theta)
  Q_emax <- comp_Q_emax(data, weight,theta)
  Q_led50 <- comp_Q_led50(data, weight,theta)
  Q <- list(Q_e0,Q_emax,Q_led50)

  K_e0 <- comp_K_e0(data, weight,theta)
  K_emax <- comp_K_emax(data, weight,theta)
  K_led50 <- comp_K_led50(data, weight,theta)
  Ka <- list(K_e0,K_emax,K_led50)

  I <- Comp_I_score(data, weight,theta)
  k <- try(solve(-I))

  if (inherits(k,'try-error')){b=c(0,0,0)}
  else{
    b <- c(0,0,0)
    for(s in 1:3){
      for(i in 1:3){
        for (j in 1:3){
          for (q in 1:3){
            b[s] <- b[s]+k[s,i]*k[j,q]*(1/2*Ka[[q]][i,j]+Q[[q]][i,j])

          }
        }
      }
    }
  }



  return(b)
}


#' Cox–Snell bias-corrected estimator (one-step using \pkg{clinDR} MLE)
#'
#' Starts from the MLE from \pkg{clinDR} \code{fitEmax} and applies a user-supplied
#' Cox–Snell bias correction.
#'
#' @param data A data.frame (or list) with \code{y} (0/1) and \code{dose}.
#' @param weight Numeric vector of case weights.
#' @param theta Numeric(3) initial guess; ignored if \pkg{clinDR} MLE succeeds.
#'
#' @details
#' Requires helpers \code{comp_bias()} and \code{Comp_I_score()} in your package.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{bias-corrected parameter vector \code{c(e0, emax, led50)}.}
#'   \item{vc}{Variance-covariance matrix (generalized inverse of information).}
#' }
#' @examples
#' \dontrun{
#' if(interactive()){
#'  theta_true=matrix(c(qlogis(0.1),qlogis(0.8)-qlogis(0.1),log(7.5)),1,3)
#'  colnames(theta_true)<- c('e_0','emax','led_50')
#'  theta_true <- as.data.frame(theta_true)
#'  dose_set <- c(0,7.5,22.5,75,225)
#'  n=355
#'  data <-sim_data(theta_true,n,dose_set)
#'  res <- comp_theta_cox_snell(data=data )
#'  }
#' }
#' @seealso
#'  \code{\link[clinDR]{fitEmax}}
#' @export
comp_theta_cox_snell<-function(data=NULL, weight=NULL, theta=NULL){
  theta_0 <-  clinDR::startEmax(y=data$y,dose=data$dose,binary = TRUE)
  theta_0 <- theta_0[c(3,2,1)]

  if(is.null(theta)){theta<-theta_0}

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  if(is.null(weight)){ weight <-rep(1,nrow(data))}
  vcv <- c()
  theta_mle <- clinDR::fitEmax(y=data$y,dose=data$dose,modType=3,binary = TRUE,diagnostics=FALSE)

  if(is.null(theta_mle)){theta_new <-c(NA,NA,NA)
  vcv <-NA}
  else{
  theta <- theta_mle$fit$estimate[c(3,2,1)]
  bias <- comp_bias(data, weight, theta=theta)
  theta_new <- theta-bias
  vcv <- try(MASS::ginv(Comp_I_score(data=data, weight=weight, theta=theta_new)))
  if (inherits(vcv,'try-error')){vcv=matrix(NA,3,3)}
}
  return(list(par=as.numeric(theta_new),vc=vcv))

}
