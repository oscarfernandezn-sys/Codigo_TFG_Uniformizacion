#' @title Calcula N usando la función de cuantil de Poisson (la forma más robusta).
#'
#' @param bt El producto de la tasa de uniformización y el tiempo (B * t).
#' @param epsilon La tolerancia de error deseada.
#'
#' @return Un número entero que representa el número de términos N requeridos.
#'
calcular_N_qpois <- function(bt, epsilon) {
  if (bt < 0) stop("bt no puede ser negativo")
  if (epsilon <= 0 || epsilon >= 1) stop("epsilon debe estar entre 0 y 1")
  
  # qpois encuentra el valor N tal que P(X <= N) >= 1 - epsilon
  # Es la forma más directa, rápida y numéricamente estable de resolver el problema.
  return(qpois(1 - epsilon, lambda = bt))
}

print(calcular_N_qpois(1, 10^-5))  
print(calcular_N_qpois(10, 10^-5)) 
print(calcular_N_qpois(100, 10^-5))  
print(calcular_N_qpois(1000, 10^-5)) 
print(calcular_N_qpois(5000, 10^-5)) 
print(calcular_N_qpois(10000, 10^-5)) 