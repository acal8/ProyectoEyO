# solucion trivial ejemplo
remove(list=ls())
#setwd("/home/emiliano/Documents/estadistica/estadistica_y_optimizacion_master/proyecto/")
dir()
# cargamos funciones con soluciones triviales. 
# Cada equipo debera definir en su archivo de solucion_teamName.R todas las
# funciones incluidas abajo

#############################################################################################
# MODIFICAR DESDE AQUI...


teamName <- "Pichones"
# integrante 1: Adrián Carrasco Alcalá
# integrante 2: Javier Herrero Pérez
# integrante 3: Clara Montalvá Barcenilla

##########################
# seccion 2 - predicciones
##########################
library(readr)
library(ggplot2)
library(forecast)
library(quadprog)
library(nloptr)

# Funcion para predecir UN activo para el periodo test en el momento t
getPred <- function(x_train, x_test_past){
  # INPUTS: x_train serie de tiempo de todo periodo train para UN activo
  #       : x_test_past serie de tiempo de periodo test que va desde 0 hast t-1
  
  # OUTPUTS: mu_hat: prediccion del activo para el momento t del periodo test
  #          se_hat: desviación estándar de la predicción. 
  
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
  # INSERTAR ESTIMACION DE MODELO DE SERIES DE TIEMPO PARA MOMENTO T PARA UN ACTIVO
  mu_hat <- # INSERTAR CALCULO FINAL DE LA PREDICCION (UN SOLO NUMERO)
  se_hat <- # INSERTAR CALCULO FINAL DE LA DESVIACIONE STANDAR (UN SOLO NUMERO)
  
  return(list(mu_hat=mu_hat, se_hat=se_hat))
}

#####################################
# seccion 3 - utilidad media-varianza
###########################################################
# 3.1 utilidad media-varianza, alfa_i positiva o negativa
###########################################################
gammaMV <- 5#INSERTAR VALOR EN REALES 

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaMV <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  R<-cor(Xpast)
  D<-diag(sig)
  Sigma<-D%*%R%*%D
  return(Sigma)
}
getSigmaDiag<-getSigmaMV
# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad Media-Varianza CON posiciones cortas
getAlphaMV <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-MV sujeto a que alpha sume a 1.
  
  #INSERTAR PASOS DE OPTIMIZACION
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

############################################
# 3.2 utilidad media-varianza, alfa_i positiva 
############################################

gammaMVPos <- 20#INSERTAR VALOR EN REALES

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaMVPos <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  
  #INSERTAR CONSTRUCCION DE MATRIZ DE COVARIANZAS 
  R<-cor(Xpast)
  D<-diag(sig)
  Sigma<-D%*%R%*%D
  return(Sigma)
}

# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad Media-Varianza SIN posiciones cortas
getAlphaMVPos <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-MV sujeto a que alpha sume a 1 y alpha_i>0.
  
  #INSERTAR PASOS DE OPTIMIZACION
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
  
  to_zero <- abs(alpha) < 1e-14
  alpha[to_zero] <- 0.0
  return(alpha)
}


################################################
# seccion 4 - 
# utilidad log, alfa_i positiva o negativa
################################################
gammaLog <- 20# INSERTAR VALOR EN REALES

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaLog <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  
  #INSERTAR CONSTRUCCION DE MATRIZ DE COVARIANZAS 
  R<-cor(Xpast)
  D<-diag(sig)
  Sigma<-D%*%R%*%D
  return(Sigma)
}

# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad log CON posiciones cortas
getAlphaLog <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-log sujeto a que alpha sume a 1.
  
  #INSERTAR PASOS DE OPTIMIZACION
  n <- length(mu)
  
  # 1. FUNCIÓN OBJETIVO (MINIMIZAR EL NEGATIVO DE LA UTILIDAD LOGARÍTMICA)
  # Usamos la sustitución de variables: solo optimizamos los primeros N-1 pesos.
  Ulog_minimizacion <- function(alpha_n_1) {
    
    # a. Calculamos el N-ésimo peso para cumplir la restricción sum(alpha) = 1
    alpha_full <- numeric(n)
    alpha_full[1:(n-1)] <- alpha_n_1
    alpha_full[n] <- 1 - sum(alpha_n_1) 
    
    # b. Calculamos el rendimiento esperado (Er)
    Er <- (t(alpha_full)%*%mu)[1,1]
    if (1 + Er <= 1e-8) { 
      return(1e10) # Devolvemos un valor muy alto para minimizar
    }    
    
    riskTerm <- (t(alpha_full) %*% Sigma %*% alpha_full) / (1 + Er)^2
    logTerm <- log(1 + Er)
    
    # Devolvemos el negativo de la utilidad (para la minimización)
    return(gamma * 0.5 * riskTerm - logTerm)
  }
  
  
  par_inicial <- rep(1/n, n-1) 
  
  sol <- optim(
    par = par_inicial,
    fn = Ulog_minimizacion,
    method = "BFGS" 
  )
  
  
  alpha_resultado <- numeric(n)
  alpha_resultado[1:(n-1)] <- sol$par
  
  # Calculamos el último peso con la restricción de suma=1
  alpha_resultado[n] <- 1 - sum(alpha_resultado[1:(n-1)])
  
  to_zero <- abs(alpha_resultado) < 1e-15 
  alpha_resultado[to_zero] <- 0.0
  
  return(alpha_resultado)
}


# ... HASTA AQUI
#############################################################################################


###############################################################
# Evaluación de soluciones
###############################################################
source("funciones/eval_funcs.R")

#setwd("/home/emiliano/Documents/estadistica/estadistica_y_optimizacion_master/proyecto/")
X <- read.csv("stock_returns_train_2.csv")
X <- ts(X/100)

# Validation mode - para que se evaluen asi mismos con el 
Xtrain <- window(X, start=1,end=8*12) # el start-end es un ejemplo, pueden cambiarlo
Xtest <- window(X, start=8*12+1,end=10*12)

# Test mode - no tendran el archivo stock_returns_test.csv asi que esto lo 
# ejecutaremos una vez entreguen soluciones
#Xtrain <- X
#Xtest <- ts(read.csv("stock_returns_test.csv"))


#seccion 2 - predicciones
set.seed(43)
res <- getPred_ts(Xtrain, Xtest, getPred)
mu_hat = res$mu_hat
se_hat = res$se_hat

# MAPE
i <- 2
plot(as.data.frame(Xtest)[,i], ty="l")
lines(mu_hat[,i], col="blue", ty="l", xlab = 'Tiempo (meses)', ylab = '')



rmse <- sqrt(mean((Xtest-mu_hat)^2))
evals <- c(rmse=rmse)
evals

# seccion 3 - utilidad media varianza
# utilidad media-varianza, alfa_i positiva o negativa

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaMV, getSigmaMV, getAlphaMV, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode="sum1")
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retMV=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaMV, getSigmaMV, passChecks, Umv)
evals <- c(evals,  Umv=Umv_rel)

# utilidad media-varianza, alfa_i positiva 

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaMVPos, getSigmaMVPos, getAlphaMVPos, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode=c("sum1","pos"))
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retMVPos=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaMVPos, getSigmaMVPos, passChecks, Umv)
evals <- c(evals,  UmvPos=Umv_rel)


# seccion 4 - 
# utilidad log, alfa_i positiva o negativa

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaLog, getSigmaLog, getAlphaLog, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode=c("sum1"))
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retLog=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaLog, getSigmaLog, passChecks, Umv)
evals <- c(evals,  UmvPosInt=Umv_rel)

evals


