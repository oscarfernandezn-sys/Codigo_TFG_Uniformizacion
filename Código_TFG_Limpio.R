# --- Cargar librerías (si son necesarias para visualización, etc.) ---
library(ggplot2) # Opcional, para gráficos
library(dplyr)

# --- FUNCIÓN 1: Determinar el punto de truncamiento M ---
# Calcula M basándose en la distribución estacionaria de la cola M/M/k.
determinar_M <- function(lambda, mu, k, epsilon1) {
  
  rho <- lambda / (k * mu)
  
  # Condición de estabilidad
  if (rho >= 1) {
    stop("El sistema es inestable (rho >= 1). 
         No se puede calcular M basado en el estado estacionario.")
  }
  
  eta <- lambda / mu
  
  # Calcular pi_0
  suma_pi0_inv <- 0
  for (n in 0:(k - 1)) {
    suma_pi0_inv <- suma_pi0_inv + (eta^n) / factorial(n)
  }
  termino_k <- (eta^k) / (factorial(k) * (1 - rho))
  pi0 <- 1 / (suma_pi0_inv + termino_k)
  
  # Calcular M iterativamente
  M <- 0
  suma_acumulada <- pi0
  pi_n <- pi0
  
  while (suma_acumulada < (1 - epsilon1)) {
    M <- M + 1
    if (M < k) {
      pi_n <- pi0 * (eta^M) / factorial(M)
    } else {
      # Para n >= k, pi_n = pi_{n-1} * rho
      # Esta es una forma más estable numéricamente que usar la fórmula 
      # directa
      pi_n <- pi_n * rho 
    }
    suma_acumulada <- suma_acumulada + pi_n
  }
  
  return(M)
}


# --- FUNCIÓN 2: Generar la matriz generadora truncada Q_M ---

generar_Q_truncada <- function(lambda, mu, k, M) {
  
  num_estados <- M + 2
  Q_M <- matrix(0, nrow = num_estados, ncol = num_estados)
  
  estados <- c(0:M, "s_a")
  rownames(Q_M) <- estados
  colnames(Q_M) <- estados
  
  for (i in 0:M) {
    i_r <- i + 1
    tasa_nacimiento <- lambda
    
    if (i > 0) {
      tasa_muerte <- min(i, k) * mu
    } else {
      tasa_muerte <- 0
    }
    
    if (i < M) {
      Q_M[i_r, i_r + 1] <- tasa_nacimiento
    } else {
      Q_M[i_r, num_estados] <- tasa_nacimiento
    }
    
    if (i > 0) {
      Q_M[i_r, i_r - 1] <- tasa_muerte
    }
    
    Q_M[i_r, i_r] <- -(tasa_nacimiento + tasa_muerte)
  }
  
  return(Q_M)
}


# --- FUNCIÓN 3: Aplicar el algoritmo de uniformización ---
# Calcula la distribución transitoria pi(t).
algoritmo_uniformizacion <- function(Q_M, pi_inicial, t, epsilon2) {
  
  # Parámetros del algoritmo
  lambda_u <- max(abs(diag(Q_M)))
  P <- diag(nrow(Q_M)) + Q_M / lambda_u
  
  # Determinar N (truncamiento de la serie)
  N <- 0
  prob_poisson <- exp(-lambda_u * t)
  suma_poisson <- prob_poisson
  
  while (suma_poisson < (1 - epsilon2)) {
    N <- N + 1
    prob_poisson <- prob_poisson * (lambda_u * t) / N
    suma_poisson <- suma_poisson + prob_poisson
  }
  
  # Calcular la suma de la serie
  suma_vectorial <- pi_inicial
  termino_vectorial <- pi_inicial
  
  for (i in 1:N) {
    termino_vectorial <- termino_vectorial %*% P
    suma_vectorial <- suma_vectorial + termino_vectorial
  }
  
  # Calcular pi(t)
  # Necesitamos un vector con las probabilidades de Poisson
  prob_poisson_vec <- dpois(0:N, lambda = lambda_u * t)
  
  # Re-cálculo más estable de pi(t)
  pi_t <- matrix(0, nrow = 1, ncol = length(pi_inicial))
  termino_vectorial <- pi_inicial
  for (i in 0:N) {
    pi_t <- pi_t + prob_poisson_vec[i+1] * termino_vectorial
    termino_vectorial <- termino_vectorial %*% P
  }
  
  # Añadir nombres de los estados al resultado
  colnames(pi_t) <- rownames(Q_M)
  
  return(list(pi_t = pi_t, N = N))
}


# --- FUNCIÓN 4: Uniformización para la matriz P(t) (Sin cambios) ---
uniformizacion_matriz <- function(Q_M, t, epsilon2) {
  lambda_u <- max(abs(diag(Q_M)))
  P <- diag(nrow(Q_M)) + Q_M / lambda_u
  
  N <- 0
  prob_poisson_term <- exp(-lambda_u * t)
  suma_poisson <- prob_poisson_term
  
  while (suma_poisson < (1 - epsilon2)) {
    N <- N + 1
    prob_poisson_term <- prob_poisson_term * (lambda_u * t) / N
    suma_poisson <- suma_poisson + prob_poisson_term
  }
  
  num_estados <- nrow(Q_M)
  P_t <- matrix(0, nrow = num_estados, ncol = num_estados)
  P_n <- diag(num_estados)
  
  for (n in 0:N) {
    prob_n <- dpois(n, lambda = lambda_u * t)
    P_t <- P_t + prob_n * P_n
    P_n <- P_n %*% P
  }
  
  return(list(P_t = P_t, N = N))
}

# --- FUNCIÓN 5: Algoritmo de Escalado y Cuadratura ---
# Realiza el cálculo sin imprimir los pasos intermedios.
escalado_y_cuadratura <- function(Q_M, pi_inicial, t, epsilon2, t0) {
  
  # 1. Escalado
  m <- ceiling(log2(t / t0))
  t_final <- t0 * 2^m
  
  # 2. Cálculo base
  resultado_base <- uniformizacion_matriz(Q_M, t0, epsilon2)
  P_t0 <- resultado_base$P_t
  
  # 3. Cuadratura
  P_final <- P_t0
  if (m > 0) {
    for (i in 1:m) {
      P_final <- P_final %*% P_final
    }
  }
  
  # 4. Resultado
  pi_t <- pi_inicial %*% P_final
  colnames(pi_t) <- rownames(Q_M)
  
  return(list(pi_t = pi_t, t_final = t_final, m = m, 
              N_base = resultado_base$N))
}


# --- FUNCIÓN 6: Calcular tiempo de convergencia ---
calcular_tiempo_convergencia <- function(lambda, mu, k, epsilon1, 
                                         epsilon2, epsilon_conv, t_max = 200) {
  # 1. Preparar el sistema truncado
  M <- determinar_M(lambda, mu, k, epsilon1)
  Q_M <- generar_Q_truncada(lambda, mu, k, M)
  pi_inicial <- c(1, rep(0, M + 1))
  
  # 2. Inicializar variables para el bucle while
  t_actual <- 0
  convergencia_alcanzada <- FALSE
  pi_t_anterior <- pi_inicial
  while (!convergencia_alcanzada && t_actual < t_max) {
    t_actual <- t_actual + 1
    umbral_t <- 20 
    if (t_actual < umbral_t) {
      resultado <- algoritmo_uniformizacion(Q_M, pi_inicial, t_actual, epsilon2)
      pi_t_actual <- resultado$pi_t
    } else {
      # Usamos un t0 dinámico basado en un divisor fijo (16 por defecto)
      divisor_fijo <- 16 
      t0_dinamico <- t_actual / divisor_fijo
      resultado <- escalado_y_cuadratura(Q_M, pi_inicial, t_actual, epsilon2, t0 = t0_dinamico)
      pi_t_actual <- resultado$pi_t
    }
    # Calcular la norma infinito de la diferencia entre la distribución actual y la anterior
    norma_inf <- max(abs(pi_t_actual - pi_t_anterior))
    # Comprobar el criterio de convergencia
    if (norma_inf < epsilon_conv) {
      convergencia_alcanzada <- TRUE 
    }
    # Actualizar para la siguiente iteración
    pi_t_anterior <- pi_t_actual
  }
  # 3. Devolver el resultado final
  if (convergencia_alcanzada) {
    cat(sprintf("Convergencia alcanzada en t = %d para caso (lambda=%.2f, mu=%.2f)\n", t_actual, lambda, mu))
    return(t_actual)
  } else {
    cat(sprintf("ADVERTENCIA: No se alcanzó convergencia en t_max = %d para caso (lambda=%.2f, mu=%.2f)\n", t_max, lambda, mu))
    return(NA)
  }
}




# --- FUNCIÓN 7: Distribucion limite ---
calcular_distribucion_limite_truncada <- function(lambda, mu, k, epsilon1) {
  
  # --- 1. Determinar dinámicamente el punto de truncamiento M ---
  M <- determinar_M(lambda, mu, k, epsilon1)
  
  # --- 2. Calcular pi_0 ---
  eta <- lambda / mu
  rho <- eta / k
  
  suma_terminos <- 0
  for (n in 0:(k - 1)) {
    suma_terminos <- suma_terminos + (eta^n) / factorial(n)
  }
  termino_k <- (eta^k / factorial(k)) * (1 / (1 - rho))
  pi_0 <- 1 / (suma_terminos + termino_k)
  
  # --- 3. Calcular pi_n para los estados 0 a M ---
  pi_n_vector <- numeric(M + 1)
  pi_n_vector[1] <- pi_0 # El índice 1 corresponde a n=0
  
  for (n in 1:M) {
    if (n < k) {
      pi_n <- (eta^n / factorial(n)) * pi_0
    } else {
      pi_n <- ((eta^k / factorial(k)) * rho^(n - k)) * pi_0
    }
    pi_n_vector[n + 1] <- pi_n
  }
  
  names(pi_n_vector) <- 0:M
  return(pi_n_vector)
}


# --- 1. Funciones Auxiliares para las Métricas ---

#' Calcula la probabilidad de que un cliente tenga que esperar (P_w(t)).

calcular_Pw <- function(pi_t, k, M) {
  # Los estados donde se espera son n = k, k+1, ..., M y el estado absorbente s_a.
  # En nuestro vector, estos corresponden a los índices desde (k+1) hasta el final.
  indices_espera <- (k + 1):length(pi_t)
  prob_espera <- sum(pi_t[indices_espera])
  return(prob_espera)
}


#' Calcula la utilización del sistema (U(t)).

calcular_U <- function(pi_t, k, M) {
  # Vector con los estados numéricos (0, 1, ..., M)
  estados_numeric <- 0:M
  # Vector que indica cuántos servidores están ocupados en cada estado 'n'
  servidores_por_estado <- pmin(estados_numeric, k)
  
  # Número esperado de servidores ocupados para los estados 0 a M
  esp_serv_0_M <- sum(servidores_por_estado * pi_t[1:(M + 1)])
  
  # Para el estado absorbente s_a (que representa n > M), todos los 'k' servidores están ocupados.
  esp_serv_sa <- k * pi_t[length(pi_t)]
  
  # Sumamos ambas partes para obtener el total esperado
  total_esp_servidores_ocupados <- esp_serv_0_M + esp_serv_sa
  
  # La utilización es el total esperado dividido entre el número de servidores
  utilizacion <- total_esp_servidores_ocupados / k
  return(utilizacion)
}

#' Calcula la probabilidad de baja utilización (<50%) (P_baja_util).
calcular_Pbaja_util <- function(pi_t, k) {
  # Límite superior de estados a sumar (n < k/2)
  # ceiling(k / 2) - 1 nos da el estado entero más alto que cumple la condición.
  limite_n <- ceiling(k / 2) - 1
  
  # Si k=1, el límite es 0. Si k=0, es -1, por lo que no se suma nada (correcto).
  if (limite_n < 0) { return(0) }
  
  # Los estados a sumar son 0, 1, ..., limite_n.
  # Sus índices en el vector pi_t son 1, 2, ..., limite_n + 1.
  indices_sumar <- 1:(limite_n + 1)
  prob_baja_util <- sum(pi_t[indices_sumar])
  return(prob_baja_util)
}

















###########################################################
# Analisis pi(t) CASO BASE CON GRAFICO
###########################################################


# --- PASO 1: CALCULAR LA DISTRIBUCIÓN LÍMITE USANDO TU FUNCIÓN ---

# Parámetros del escenario base (asegúrate de que estén definidos)
lambda <- 8
mu <- 1
k <- 9
epsilon1 <- 1e-5 # Tolerancia para determinar M
M <- determinar_M(lambda, mu, k, epsilon1)
Q_M <- generar_Q_truncada(lambda, mu, k, M)
pi_inicial <- c(1, rep(0, M + 1))

# Llamamos a tu función para obtener el vector con la distribución límite
pi_limite_vector <- calcular_distribucion_limite_truncada(lambda, mu, k, epsilon1)

# Convertimos el vector resultado en un data.frame manejable para ggplot
datos_limite <- data.frame(
  Estado = as.numeric(names(pi_limite_vector)),
  Prob_Limite = as.vector(pi_limite_vector)
)


# --- PASO 2: CALCULAR LA DISTRIBUCIÓN TRANSITORIA Y UNIR LOS DATOS ---

# Tiempos para los paneles del gráfico
tiempos_paneles <- c(2, 8, 16, 64)
epsilon2 <- 1e-5
lista_de_resultados_pi <- list()


for (i in 1:length(tiempos_paneles)) {
  t_panel <- tiempos_paneles[i]
  if (t_panel >= 30) {
    resultado <- escalado_y_cuadratura(Q_M, pi_inicial, t_panel, epsilon2, t0 = t_panel / 16)
  } else {
    resultado <- algoritmo_uniformizacion(Q_M, pi_inicial, t_panel, epsilon2)
  }
  
  lista_de_resultados_pi[[i]] <- data.frame(
    t = t_panel,
    Estado = 0:M,
    Probabilidad = as.vector(resultado$pi_t)[1:(M + 1)]
  )
}

datos_grafico_pi_transitoria <- do.call(rbind, lista_de_resultados_pi)

# Unimos los dos dataframes usando la columna "Estado" como clave
datos_grafico_pi <- merge(datos_grafico_pi_transitoria, datos_limite, by = "Estado")

# Preparamos los datos finales para el gráfico
datos_grafico_pi$Panel_Titulo <- factor(paste("t =", datos_grafico_pi$t, "minutos"), levels = paste("t =", tiempos_paneles, "minutos"))
max_prob <- max(c(datos_grafico_pi$Probabilidad, datos_grafico_pi$Prob_Limite))
y_lim_max <- ceiling(max_prob * 20) / 20


# --- PASO 3: CÓDIGO DEL GRÁFICO (AHORA FUNCIONA CORRECTAMENTE) ---

grafico_cuadricula_pi <- ggplot(data = datos_grafico_pi, aes(x = Estado)) +
  geom_col(aes(y = Probabilidad, fill = "transitoria"), color = "black", width = 0.9) +
  geom_line(aes(y = Prob_Limite, color = "limite"), linewidth = 0.6) +
  
  scale_fill_manual(
    name = NULL,
    values = c("transitoria" = "black"),
    labels = expression(paste("Distribución transitoria ", pi(t)))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("limite" = "red"),
    labels = expression(paste("Distribución límite teórica ", pi))
  ) +
  
  facet_wrap(~ Panel_Titulo, nrow = 1) +
  coord_cartesian(xlim = c(-0.5, 100), ylim = c(0, y_lim_max)) +
  labs(
    title = "Evolución de la Distribución de Probabilidad π(t)",
    subtitle = "El sistema comienza en el estado 0",
    x = "Número de Clientes en el Sistema (n)",
    y = "Probabilidad P(X(t)=n)",
    fill = "",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )

# Visualizar el gráfico
print(grafico_cuadricula_pi)




###########################################################
# Analisis sensibilidad
###########################################################

# --- 1. Definición de Todos los Escenarios de Análisis---

plan_analisis <- data.frame(
  caso_id = c("A1", "A2", "A3 (Base)", "A4", "A5", 
              "B1", "B2", "B3 (Base)", "B4", "B5", 
              "C1", "C2", "C3 (Base)", "C4", "C5"),
  escenario_grupo = c(rep("A. Tasa de Llegada (β)", 5), 
                      rep("B. Tasa de Servicio (μ)", 5), 
                      rep("C. Nº de Servidores (k)", 5)),
  lambda = c(7.1, 7.55, 8.0, 8.45, 8.8, 
             8.0, 8.0, 8.0, 8.0, 8.0, 
             8.0, 8.0, 8.0, 8.0, 8.0),
  mu =   c(1.0, 1.0, 1.0, 1.0, 1.0,
           0.92, 0.95, 1.0, 1.05, 1.1,
           1.0, 1.0, 1.0, 1.0, 1.0),
  k =    c(9, 9, 9, 9, 9,
           9, 9, 9, 9, 9,
           7, 8, 9, 10, 11)
)


# --- Bucle de Cálculo de Sensibilidad ---
# Tarda aproximadamente 7 minutos desde que comienza hasta que genera
# todos los gráficos
tiempos_analisis_sens <- seq(0, 64, by = 1)
resultados_completos_df <- data.frame()
cat("\nINICIANDO ANÁLISIS DE SENSIBILIDAD COMPLETO...\n")
for (i in 1:nrow(plan_analisis)) {
  caso_actual <- plan_analisis[i, ]
  cat(paste0("\n--- PROCESANDO CASO: ", caso_actual$caso_id, " ---\n"))
  tryCatch({
    M_sens <- determinar_M(caso_actual$lambda, caso_actual$mu, caso_actual$k, 1e-5)
    Q_M_sens <- generar_Q_truncada(caso_actual$lambda, caso_actual$mu, caso_actual$k, M_sens)
    pi_inicial_sens <- c(1, rep(0, M_sens + 1)) 
    resultados_caso_df <- data.frame()
    for (t_actual in tiempos_analisis_sens) {
      if (t_actual == 0) {
        pi_t_actual <- pi_inicial_sens
        t_final_real <- 0
      } else if (t_actual < 30) {
        resultado <- algoritmo_uniformizacion(Q_M_sens, pi_inicial_sens, t_actual, 1e-5)
        pi_t_actual <- resultado$pi_t
        t_final_real <- t_actual
      } else { # t >= 30
        resultado <- escalado_y_cuadratura(Q_M_sens, pi_inicial_sens, t_actual, 1e-5, t0 = t_actual / 16)
        pi_t_actual <- resultado$pi_t
        t_final_real <- resultado$t_final
      }
      resultados_caso_df <- rbind(resultados_caso_df, data.frame(
        t = t_final_real,
        Pw = calcular_Pw(pi_t_actual, caso_actual$k, M_sens),
        U = calcular_U(pi_t_actual, caso_actual$k, M_sens),
        Pbaja_util = calcular_Pbaja_util(pi_t_actual, caso_actual$k)
      ))
    }
    resultados_caso_df$caso_id <- caso_actual$caso_id
    resultados_completos_df <- rbind(resultados_completos_df, resultados_caso_df)
  }, error = function(e) {
    cat(paste("ERROR en el caso", caso_actual$caso_id, ":", e$message, "-> Se omite este caso.\n"))
  })
}
cat("\n...ANÁLISIS DE SENSIBILIDAD FINALIZADO.\n\n")

# --- Gráficos de Sensibilidad ---
# --- Función de Gráficos Comparativos (con colores neutros) ---

crear_grafico_comparativo <- function(datos, casos_a_filtrar, metrica, titulo, subtitulo, eje_y) {
  
  datos_filtrados <- datos %>% 
    filter(caso_id %in% casos_a_filtrar) %>%
    # Aseguramos el orden correcto de las leyendas
    mutate(caso_id = factor(caso_id, levels = casos_a_filtrar)) 
  
  # Paleta de colores neutros. Puedes cambiar estos códigos hexadecimales si lo deseas.
  paleta_neutra <- c("#000000", "#1E6F5C", "#2940D3", "#CC0000", "#607EAA")
  
  grafico <- ggplot(data = datos_filtrados, aes(x = t, y = .data[[metrica]], color = caso_id, group = caso_id)) +
    geom_line(linewidth = 1) +
    
    # LÍNEA AÑADIDA: Asigna manualmente los colores de nuestra paleta
    scale_color_manual(values = paleta_neutra) +
    
    labs(
      title = titulo,
      subtitle = subtitulo,
      x = "Tiempo (minutos)",
      y = eje_y,
      color = "Caso" # Título de la leyenda
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      # --- LÍNEAS AÑADIDAS PARA LA LEYENDA ---
      legend.title = element_text(size = 16, face = "bold"), # Aumenta y pone en negrita el título "Caso"
      legend.text = element_text(size = 16)                 # Aumenta el tamaño del texto de las etiquetas
    )
  
  
  return(grafico)
}

# --- 3. Generación de Gráficos Comparativos (Versión 5x3) ---

# (La función crear_grafico_comparativo no necesita cambios)

# --- GRÁFICOS PARA ESCENARIO A (variando lambda) ---
casos_A <- c("A1", "A2", "A3 (Base)", "A4", "A5")
grafico_A_Pw <- crear_grafico_comparativo(resultados_completos_df, casos_A, "Pw", "Sensibilidad de Pw a la Tasa de Llegada (λ)", "μ=1, k=9", "Prob. de Espera (Pw)")
grafico_A_U <- crear_grafico_comparativo(resultados_completos_df, casos_A, "U", "Sensibilidad de U a la Tasa de Llegada (λ)", "μ=1, k=9", "Utilización (U)")
grafico_A_Pbaja <- crear_grafico_comparativo(resultados_completos_df, casos_A, "Pbaja_util", "Sensibilidad de P_baja_util a la Tasa de Llegada (β)", "μ=1, k=9", "Prob. Baja Utilización")

# --- GRÁFICOS PARA ESCENARIO B (variando mu) ---
casos_B <- c("B1", "B2", "B3 (Base)", "B4", "B5")
grafico_B_Pw <- crear_grafico_comparativo(resultados_completos_df, casos_B, "Pw", "Sensibilidad de Pw a la Tasa de Servicio (μ)", "λ=8, k=9", "Prob. de Espera (Pw)")
grafico_B_U <- crear_grafico_comparativo(resultados_completos_df, casos_B, "U", "Sensibilidad de U a la Tasa de Servicio (μ)", "λ=8, k=9", "Utilización (U)")
grafico_B_Pbaja <- crear_grafico_comparativo(resultados_completos_df, casos_B, "Pbaja_util", "Sensibilidad de P_baja_util a la Tasa de Servicio (μ)", "β=8, k=9", "Prob. Baja Utilización")

# --- GRÁFICOS PARA ESCENARIO C (variando k) ---
# Omitimos C1 y C2 porque son inestables y no generan datos
casos_C <- c("C3 (Base)", "C4", "C5") 
grafico_C_Pw <- crear_grafico_comparativo(resultados_completos_df, casos_C, "Pw", "Sensibilidad de Pw al Número de Servidores (k)", "λ=8, μ=1", "Prob. de Espera (Pw)")
grafico_C_U <- crear_grafico_comparativo(resultados_completos_df, casos_C, "U", "Sensibilidad de U al Número de Servidores (k)", "λ=8, μ=1", "Utilización (U)")
grafico_C_Pbaja <- crear_grafico_comparativo(resultados_completos_df, casos_C, "Pbaja_util", "Sensibilidad de P_baja_util al Número de Servidores (k)", "β=8, μ=1", "Prob. Baja Utilización")

# Mostramos todos los gráficos generados
print(grafico_A_Pw); print(grafico_A_U); print(grafico_A_Pbaja)
print(grafico_B_Pw); print(grafico_B_U); print(grafico_B_Pbaja)
print(grafico_C_Pw); print(grafico_C_U); print(grafico_C_Pbaja)




###########################################################
# Análisis convergencia
###########################################################
# ESTA PRIMERA PARTE HASTA EL ANALISIS AMPLIADO (SIN CONTARLO)
# TARDA 12 MINUTOS

# CASOS A PARA EL CUADRO 3.2
# Filtrar para seleccionar ÚNICAMENTE los escenarios del grupo A
escenarios_convergencia <- plan_analisis %>% 
  filter(startsWith(caso_id, "A"))

# Crear data frame para almacenar los resultados
resultados_convergencia_df <- data.frame()

for (i in 1:nrow(escenarios_convergencia)) {
  caso_actual <- escenarios_convergencia[i, ]
  
  cat(paste0("--- Calculando convergencia para CASO: ", caso_actual$caso_id, " (λ=", caso_actual$lambda, ", μ=", caso_actual$mu, ") ---\n"))
  
  tryCatch({
    tiempo_conv <- calcular_tiempo_convergencia(
      lambda = caso_actual$lambda,
      mu = caso_actual$mu,
      k = caso_actual$k,
      epsilon1 = 1e-5,
      epsilon2 = 1e-5,
      epsilon_conv = 1e-4/2
    )
    
    resultados_convergencia_df <- rbind(resultados_convergencia_df, data.frame(
      caso_id = caso_actual$caso_id,
      lambda = caso_actual$lambda,
      rho = round(caso_actual$lambda / (caso_actual$k * caso_actual$mu), digits = 2),
      tiempo_convergencia = tiempo_conv
    ))
    
  }, error = function(e) {
    cat(paste("ERROR en el caso", caso_actual$caso_id, ":", e$message, "\n"))
  })
}

cat("\n...ANÁLISIS DE CONVERGENCIA FINALIZADO.\n\n")

# --- Presentación de Resultados ---

# 1. Mostrar la tabla de resultados
cat("Tabla de Resultados: Tiempo de Convergencia (t*) para Escenarios A\n")
print(resultados_convergencia_df)



###########################################################
# Análisis convergencia AMPLIADO
###########################################################
# TARDA EN TORNO A 30 MINUTOS


# Definir los parámetros fijos
k_fijo <- 9
mu_fijo <- 1.0
# Crear una secuencia de 15 valores de rho equidistantes
rho_seq <- seq(from = 0.725, to = 0.985, length.out = 15)
# Calcular los valores de lambda correspondientes para cada rho
lambda_seq <- rho_seq * k_fijo * mu_fijo
# Crear el data frame con los nuevos escenarios
escenarios_convergencia <- data.frame(
  lambda = lambda_seq,
  mu = mu_fijo,
  k = k_fijo
)
escenarios_convergencia$caso_id <- paste0("A", 1:nrow(escenarios_convergencia))
# Preparar el data frame para los resultados
resultados_convergencia_df <- data.frame()
# Ejecutar el bucle para cada escenario
for (i in 1:nrow(escenarios_convergencia)) {
  caso_actual <- escenarios_convergencia[i, ]
  cat(paste0("--- Calculando convergencia para CASO: ", caso_actual$caso_id, " (lambda=", caso_actual$lambda, ", mu=", caso_actual$mu, ") ---\n"))
  tryCatch({
    tiempo_conv <- calcular_tiempo_convergencia(
      lambda = caso_actual$lambda,
      mu = caso_actual$mu,
      k = caso_actual$k,
      epsilon1 = 1e-5,
      epsilon2 = 1e-5,
      epsilon_conv = 1e-4/2
    )
    resultados_convergencia_df <- rbind(resultados_convergencia_df, data.frame(
      caso_id = caso_actual$caso_id,
      lambda = caso_actual$lambda,
      rho = rho_seq[i],
      tiempo_convergencia = tiempo_conv
    ))
  }, error = function(e) {
    cat(paste("ERROR en el caso", caso_actual$caso_id, ":", e$message, "\n"))
  })
}
cat("\n...ANÁLISIS DE CONVERGENCIA FINALIZADO.\n\n")

# --- 3. Crear el gráfico de convergencia vs. intensidad de tráfico ---
# FIGURA 3.7 EN EL TFG
# (El código del gráfico es el mismo, solo que ahora el eje X tendrá valores equidistantes)
grafico_convergencia <- ggplot(data = resultados_convergencia_df, aes(x = rho, y = tiempo_convergencia)) +
  geom_line(color = "black", linetype = "solid") +
  geom_point(color = "black", size = 1.5, alpha = 0.8) +
  labs(
    title = "Tiempo de Convergencia en funcion de la intensidad de tráfico",
    x = "Intensidad de Tráfico (ρ)",
    y = "Tiempo de Convergencia (t*)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  xlim(0.7, 1) + # Esta línea se ha añadido
  ylim(0, max(resultados_convergencia_df$tiempo_convergencia, na.rm = TRUE) * 1.05)

print(grafico_convergencia)

