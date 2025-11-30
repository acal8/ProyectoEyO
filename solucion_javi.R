library(readr)
library(ggplot2)
library(forecast)
library(quadprog)
# seccion 2 - predicciones

getPred <- function(x_train, x_test_past){
  mod <- auto.arima(x_train,stepwise =F,approximation = F)
  x_past = c(x_train, x_test_past)
  fit_upd <- Arima(x_past, model = mod)
  fc <- forecast::forecast(fit_upd, h = 1)
  mu_hat  <- as.numeric(fc$mean[1])
  z80 <- qnorm(0.8)
  se_from_up <- (fc$upper[,"80%"][1] - fc$mean[1]) / z80
  se_from_lo <- (fc$mean[1]  - fc$lower[,"80%"][1]) / z80
  se_hat_aux <- pmax(se_from_up, se_from_lo) 
  se_hat  <- as.numeric(se_hat_aux)
  return(list(mu_hat=mu_hat, se_hat=se_hat))
}

# seccion 3 - utilidad media-varianz

getSigma <- function(sig, Xpast){
  # sig es se_hat que se saca de getPred
  # Hay que calcular \Sigma_t=DRD
  # D viene definido por la diagonal de sig
  # R es una matriz de correlación de los valores pasadso Xpast
  R<-cor(Xpast)
  D<-diag(sig)
  sigma_hat<-D%*%R%*%D
  return(Sigma_hat)
}


# alpha random
getAlphaMV <- function(mu,Sigma, gamma){
  # Calculado a partir de la expresión analítica que he calculado en papel 
  n<-length(mu)
  I<-matrix(1,nrow=n,ncol=1)
  Sigma_inversa <- solve(Sigma)
  lambda<- (t(I)%*%Sigma_inversa%*%mu-gamma)/(t(I)%*%Sigma_inversa%*%I)
  alpha <- Sigma_inversa%*%(mu-lambda%*%I)/gamma
  return(alpha)
}





# alpha random y positivo
getAlphaMVPos <- function(mu,Sigma, gamma){
  # Usar la librería quadprog
  
  n<-length(mu)
  
  Dmat_qp <- gamma * Sigma
  dvec_qp <- mu
  
  
  # A. Restricción de IGUALDAD (Sum(alpha) = 1)
  A_eq <- matrix(1, nrow = n, ncol = 1)
  b_eq <- 1
  
  
  # B. Restricciones de DESIGUALDAD (alpha_i >= 0)
  A_ineq <- diag(n) # Matriz Identidad
  b_ineq <- rep(0, n)
  
  
  Amat_qp<-cbind(A_eq, A_ineq)
  bvec_qp <- c(b_eq, b_ineq)
  dval<-mu # Para seguir la expresión de la libería
  
  solucion <- solve.QP(
    Dmat = Dmat_qp,
    dvec = dvec_qp,
    Amat = Amat_qp,
    bvec = bvec_qp,
    meq = 1 # Primera restricción
  )
  alpha<- solucion$solution
  return(alpha)
}

# seccion 4 - 
# utilidad log, alfa_i positiva o negativa

getAlphaLog <- function(mu,Sigma, gamma){
  # Calculado a partir de la expresión analítica que he calculado en papel 
  n<-length(mu)
  I<-matrix(1,nrow=n,ncol=1)
  Sigma_inversa <- solve(Sigma)
  lambda<- (t(I)%*%Sigma_inversa%*%mu-gamma)/(t(I)%*%Sigma_inversa%*%I)
  alpha <- Sigma_inversa%*%(mu-lambda%*%I)/gamma
  return(alpha)
}
