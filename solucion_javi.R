library(readr)
library(ggplot2)
library(forecast)
library(quadprog)
library(nloptr)
# seccion 2 - predicciones

getPred000 <- function(x_train, x_test_past){
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

getSigmaDiag <- function(sig, Xpast){
  # sig es se_hat que se saca de getPred
  # Hay que calcular \Sigma_t=DRD
  # D viene definido por la diagonal de sig
  # R es una matriz de correlación de los valores pasadso Xpast
  R<-cor(Xpast)
  D<-diag(sig)
  Sigma_hat<-D%*%R%*%D
  return(Sigma_hat)
}


# alpha random
getAlphaRandom <- function(mu,Sigma, gamma){
  # Calculado a partir de la expresión analítica que he calculado en papel 
  n<-length(mu)
  I<-matrix(1,nrow=n,ncol=1)
  Sigma_inversa <- solve(Sigma)
  numerador <- t(I) %*% Sigma_inversa %*% mu - gamma
  denominador <- t(I) %*% Sigma_inversa %*% I
  lambda <- as.numeric(numerador / denominador)
  
  alpha <- Sigma_inversa%*%(mu-lambda*I) / gamma
  
  to_zero <- abs(alpha) < 1e-15
  alpha[to_zero] <- 0.0
  
  return(alpha)
}


# alpha random y positivo
getAlphaRandomPos <- function(mu,Sigma, gamma){
  # Usar la librería quadprog
  
  n<-length(mu)
  Dmat_qp <- gamma*Sigma
  dvec_qp <- mu
  
  # Restricción de igualdad Sum(alpha) = 1
  A_eq <- matrix(1, nrow = n, ncol = 1)
  b_eq <- 1
  
  # Restricciones de desigualdad alpha_i >= 0
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
  
  to_zero <- abs(alpha) < 1e-15
  alpha[to_zero] <- 0.0
  
  return(alpha)
}



# seccion 4 - 
# utilidad log, alfa_i positiva
# OJO! ANTES ESTABA PUESTO PARA ALFA_I POSITIVA O NEGATIVA, PERO EN EL ENUNCIADO PONE
# QUE TIENE QUE SER SIN POSICIONES CORTAS (ALPHA>=0)
getAlphaLog <- function(mu, Sigma, gamma){
  n <- length(mu)
  
  Ulog_minimizacion <- function(alpha) {
    Er <- (t(alpha)%*%mu)[1,1]
    
    if (1 + Er <= 1e-10) {
      return(1e20)
    }
    
    riskTerm <- (t(alpha)%*%Sigma%*%alpha)[1,1]/((1+Er)^2)
    logTerm <- log(1 + Er)
    
    # Multiplicamos el término de riesgo por gamma
    return(gamma*0.5*riskTerm - logTerm) 
  }
  
  # Restricción de igualdad Sum(alpha) = 1
  eval_g_eq <- function(alpha) {
    return(sum(alpha) - 1) 
  }
  
  # Restricciones de desigualdad alpha_i >= 0
  # En la función NLOPTR se necesita h(alpha) <= 0
  eval_g_ineq <- function(alpha) {
    return(-alpha)
  }
  
  # Punto de partida
  par_inicial <- rep(1/n, n) 
  
  solucion <- nloptr(
    x0 = par_inicial,
    eval_f = Ulog_minimizacion,
    eval_g_eq = eval_g_eq,
    eval_g_ineq = eval_g_ineq,
    opts = list(algorithm = "NLOPT_LN_COBYLA", 
                xtol_rel = 1e-6)
  )
  
  alpha_resultado <- solucion$solution
  
  
  to_zero <- abs(alpha_resultado) < 1e-15
  alpha_resultado[to_zero] <- 0.0
  
  return(alpha_resultado)
}


