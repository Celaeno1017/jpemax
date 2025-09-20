
comp_P_e0 <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))

    I <- I+I_i
  }

  I

}

comp_P_emax <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))*data$dose[i]/(exp(led50)+data$dose[i])

    I <- I+I_i
  }

  I

}


comp_P_led50 <- function(data, weight, theta){

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

    I_i <- I_i*weight[i]^2*(pr[i]*(1-pr[i])^3-pr[i]^3*(1-pr[i]))*(-emax*data$dose[i]*exp(led50)/(exp(led50)+data$dose[i])^2)

    I <- I+I_i
  }

  I

}





Score_e0 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q1 <- sum(weight*(data$y-pr))
  Q_e0 <- comp_Q_e0(data, weight, theta)
  P_e0 <- comp_P_e0(data, weight, theta)
  I <- Comp_I(data, weight, theta)

  I_inv <- try(MASS::ginv(I))
  if (inherits(I_inv,'try-error')){I_inv=matrix(NA,3,3)}
  U <- q1+sum(diag(I_inv%*%(Q_e0+P_e0)))/2

  U
}

Score_emax <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q2 <- sum(weight*(data$y-pr)*data$dose/(exp(led50)+data$dose))
  Q_emax <- comp_Q_emax(data, weight, theta)
  P_emax <- comp_P_emax(data, weight, theta)
  I <- Comp_I(data, weight, theta)
  I_inv <- try(MASS::ginv(I))
  if (inherits(I_inv,'try-error')){I_inv=matrix(NA,3,3)}

  U <- q2+sum(diag(I_inv%*%(Q_emax+P_emax)))/2

  U
}

Score_led50 <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  pr <- plogis(emax*data$dose/(exp(led50)+data$dose)+e0)
  q3 <- sum(weight*(data$y-pr)*(-1)*data$dose*emax*exp(led50)/(exp(led50)+data$dose)^2)
  Q_led50 <- comp_Q_led50(data, weight, theta)
  P_led50 <- comp_P_led50(data, weight, theta)
  I <- Comp_I(data, weight, theta)
  I_inv <- try(MASS::ginv(I))
  if (inherits(I_inv,'try-error')){I_inv=matrix(NA,3,3)}

  U <- q3+sum(diag(I_inv%*%(Q_led50+P_led50)))/2

  U
}

Comp_Fish_inf <- function(data, weight, theta){
  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  H <- matrix(0,3,3)

  H[1,] <- numDeriv::grad(func=Score_e0,x=theta,data=data, weight=weight)
  H[2,] <- numDeriv::grad(func = Score_emax,x=theta,data=data, weight=weight)
  H[3,] <- numDeriv::grad(func = Score_led50,x=theta,data=data, weight=weight)


  I_inv <- try(MASS::ginv(-H))
  if (inherits(I_inv,'try-error')){I_inv=matrix(NA,3,3)}
  I_inv
}





#' Firth-corrected estimating equation solution (score-based)
#'
#' Solves the Firth-corrected score equations for an Emax-type binary-response model.
#' Requires user-provided score and information helpers.
#'
#' @param data A data.frame (or list) with at least \code{y} (0/1) and \code{dose}.
#' @param weight Numeric vector of case weights, same length as \code{data$y}.
#' @param theta Numeric vector (length 3): \code{c(e0, emax, led50)} for initialization.
#'
#' @details
#' This function depends on internal helpers you must supply in the package:
#' \code{Score_e0()}, \code{Score_emax()}, \code{Score_led50()},
#' \code{Comp_Fish_inf()}, and \code{Comp_I_score()}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{par}{numeric(3) estimated parameters.}
#'   \item{Fisher.inf}{Observed/expected information matrix at the solution (from \code{Comp_I_score}).}
#'   \item{vc}{Fisher-based variance/covariance (from \code{Comp_Fish_inf}).}
#'   \item{score}{Score vector evaluated at the solution.}
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
#'  res <- comp_theta_firth_score(data=data )
#'  }
#' }
#' @seealso
#'  \code{\link[clinDR]{fitEmax}}
#' @export
comp_theta_firth_score<-function(data=NULL, weight=NULL, theta=NULL){
  if(is.null(weight)){ weight <-rep(1,nrow(data))}
  theta_0 <-  clinDR::startEmax(y=data$y,dose=data$dose,binary = TRUE)
  theta_0 <- theta_0[c(3,2,1)]

  if(is.null(theta)){theta<-theta_0}

  theta <- as.numeric(theta)
  e0 <- theta[1]
  emax <- theta[2]
  led50 <- theta[3]

  # for(i in 1:100){
  #   score <- c(Score_e0(data, weight, theta),Score_emax(data, weight, theta),
  #              Score_led50(data, weight, theta))
  #   Fish.Inv <- Comp_Fish_inf(data, weight, theta)
  #   #Fish.Inv <- solve(Comp_I(data, weight, theta))
  #
  #   theta =theta + 0.1*t(score)%*%Fish.Inv
  #
  # }
  #

  target <- function(theta,data,weight)
  {

    y <- numeric(3)
    y[1] <- Score_e0(data=data,theta = theta,weight = weight)
    y[2] <- Score_emax(data=data,theta = theta,weight = weight)
    y[3] <- Score_led50(data=data,theta = theta,weight = weight)
    y
  }

  # theta_new <-nleqslv(x=theta,fn=target,a=0,b=0,c=0,data=data,weight=weight,method = 'Newton',
  #                 jacobian = TRUE)

  theta_new <-BB::dfsane(par =theta,fn = target, control = list(trace = FALSE,maxit = 1000),data=data,weight=weight)

  I <-Comp_I_score(data=data,weight=weight,theta=theta_new$par)
  vcv <- Comp_Fish_inf(data=data,weight=weight,theta=theta_new$par)
  score <- c(Score_e0(data, weight, theta_new$par),Score_emax(data, weight, theta_new$par),
             Score_led50(data, weight, theta_new$par))
  vc <- try(MASS::ginv(I))
  if (inherits(vc,'try-error')){vc=matrix(NA,3,3)}
  return(list(par=as.numeric(theta_new$par),Fisher.inf=I, vc=vc,score=score))

}
