library(readr)
library(ggplot2)
library(forecast)
library(quadprog)
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
  return(alpha)
}


# alpha random y positivo
getAlphaRandomPos <- function(mu,Sigma, gamma){
  # Usar la librería quadprog
  
  n<-length(mu)
  
  Dmat_qp <- gamma*Sigma
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


# alpha random y positivo y en enteros/5 NO HACER!!!!
getAlphaRandomPosInt <- function(mu, Sigma, gamma){
  smpl <- c()
  sumSmpl <- 0
  for(i in 1:5){
    smpl <- c(smpl, sample(0:(5-sumSmpl),1))
    sumSmpl <- sum(smpl)
  }
  alpha <- smpl/5
  return(alpha)
}



# seccion 4 - 
# utilidad log, alfa_i positiva o negativa

getAlphaLog <- function(mu, Sigma, gamma){
  
  n <- length(mu)
  
  # --- 1. Función Objetivo a MINIMIZAR (-Ulog) ---
  # El argumento 'alpha_n_1' son los n-1 pesos que optimiza 'optim'.
  Ulog_minimizacion <- function(alpha_n_1) {
    
    # Reconstrucción del vector completo (5 variables)
    alpha_full <- numeric(n)
    alpha_full[1:(n-1)] <- alpha_n_1
    alpha_full[n] <- 1 - sum(alpha_n_1) # Restricción de suma: alpha[n] = 1 - Sum(alpha[1:4])
    
    # Calculamos Er y evitamos log(<=0)
    Er <- sum(alpha_full * mu)
    
    if (1 + Er <= 0) {
      return(1e10) # Penalización si el rendimiento esperado es demasiado negativo
    }
    
    # 2. Fórmula de -Ulog (Función a minimizar)
    # Ulog = log(1+Er) - 0.5 * (alpha^T Sigma alpha) / (1+Er)^2
    # -Ulog = 0.5 * (alpha^T Sigma alpha) / (1+Er)^2 - log(1+Er)
    
    riskTerm <- (t(alpha_full) %*% Sigma %*% alpha_full) / (1 + Er)^2
    logTerm <- log(1 + Er)
    
    # Devolvemos el negativo de la utilidad
    return(0.5 * riskTerm - logTerm)
  }
  
  # --- 2. Parámetros de Optimización ---
  par_inicial <- rep(1/n, n-1) # Punto de partida: [0.2, 0.2, 0.2, 0.2]
  
  # Llamamos a optim
  # Pasamos los datos fijos (mu, Sigma, gamma) usando '...'
  sol <- optim(
    par = par_inicial,
    fn = Ulog_minimizacion,
    method = "BFGS" # Adecuado para optimización no lineal sin restricciones de límites
  )
  
  # --- 3. Construcción del Resultado Final ---
  alpha_resultado <- numeric(n)
  
  # Asignamos los parámetros óptimos
  alpha_resultado[1:(n-1)] <- sol$par
  
  # Calculamos el último peso con la restricción de suma=1
  alpha_resultado[n] <- 1 - sum(alpha_resultado[1:(n-1)])
  
  # La normalización final 'alpha_resultado / sum(alpha_resultado)' es redundante
  # ya que el algoritmo está diseñado para que la suma sea 1.
  
  return(alpha_resultado)
}
