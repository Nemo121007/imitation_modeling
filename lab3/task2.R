# Функция для чтения данных
read_data <- function(filename) {
  data <- read.table(filename, header = FALSE)
  colnames(data) <- c("doc", "word", "count")
  return(data)
}

# Основная функция Collapsed Gibbs Sampling
collapsed_gibbs_lda <- function(data, T, alpha, beta, n_iter = 1000) {
  # Извлечение информации из данных
  D <- max(data$doc)  # Количество документов
  W <- max(data$word) # Размер словаря
  N <- nrow(data)     # Количество уникальных пар (документ, слово)
  
  # Разворачиваем данные: каждое слово повторяется count раз
  docs <- rep(data$doc, data$count)
  words <- rep(data$word, data$count)
  N_total <- sum(data$count) # Общее число слов в корпусе
  
  # Инициализация тем случайным образом
  z <- sample(1:T, N_total, replace = TRUE)
  
  # Счётчики
  NWT <- matrix(0, nrow = W, ncol = T) # Число слов w в теме t
  NTD <- matrix(0, nrow = T, ncol = D) # Число слов в теме t в документе d
  NT <- numeric(T)                     # Общее число слов в теме t
  ND <- numeric(D)                     # Общее число слов в документе d
  
  # Инициализация счётчиков
  for (i in 1:N_total) {
    d <- docs[i]
    w <- words[i]
    t <- z[i]
    NWT[w, t] <- NWT[w, t] + 1
    NTD[t, d] <- NTD[t, d] + 1
    NT[t] <- NT[t] + 1
    ND[d] <- ND[d] + 1
  }
  
  # Collapsed Gibbs Sampling
  for (iter in 1:n_iter) {
    for (i in 1:N_total) {
      d <- docs[i]
      w <- words[i]
      t <- z[i]
      
      # Уменьшаем счётчики для текущего слова
      NWT[w, t] <- NWT[w, t] - 1
      NTD[t, d] <- NTD[t, d] - 1
      NT[t] <- NT[t] - 1
      ND[d] <- ND[d] # ND не меняется, так как это общее число слов в документе
      
      # Вычисляем вероятности для каждой темы
      probs <- numeric(T)
      for (t in 1:T) {
        # Первая часть: p(z_{d,n} = t | ...)
        p1 <- (alpha + NTD[t, d]) / (alpha * T + ND[d])
        # Вторая часть: p(w_{d,n} = w | z_{d,n} = t, ...)
        p2 <- (beta + NWT[w, t]) / (beta * W + NT[t])
        probs[t] <- p1 * p2
      }
      
      # Нормализация вероятностей
      probs <- probs / sum(probs)
      
      # Сэмплируем новую тему
      t_new <- sample(1:T, 1, prob = probs)
      z[i] <- t_new
      
      # Обновляем счётчики
      NWT[w, t_new] <- NWT[w, t_new] + 1
      NTD[t_new, d] <- NTD[t_new, d] + 1
      NT[t_new] <- NT[t_new] + 1
    }
  }
  
  # Оценка theta и phi
  theta <- matrix(0, nrow = D, ncol = T)
  phi <- matrix(0, nrow = T, ncol = W)
  
  for (d in 1:D) {
    for (t in 1:T) {
      theta[d, t] <- (alpha + NTD[t, d]) / (alpha * T + ND[d])
    }
  }
  
  for (t in 1:T) {
    for (w in 1:W) {
      phi[t, w] <- (beta + NWT[w, t]) / (beta * W + NT[t])
    }
  }
  
  return(list(theta = theta, phi = phi, z = z))
}

# Тестирование на двух наборах данных
# Тест 1
data1 <- read_data("test1.dat")
result1 <- collapsed_gibbs_lda(data1, T = 3, alpha = 1, beta = 1, n_iter = 1000)
cat("Test 1:\n")
cat("Theta:\n")
print(result1$theta)
cat("Phi:\n")
print(result1$phi)
# Проверка для theta
print(rowSums(result1$theta))  # Должно быть вектор единиц длиной D
# Проверка для phi
print(rowSums(result1$phi))    # Должно быть вектор единиц длиной T

# Тест 2
data2 <- read_data("test2.dat")
result2 <- collapsed_gibbs_lda(data2, T = 20, alpha = 0.1, beta = 0.1, n_iter = 10)
cat("Test 2:\n")
cat("Theta:\n")
theta_2 <- result2$theta
phi_2 <- result2$phi
print(result2$theta)
cat("Phi:\n")
print(result2$phi)
# Проверка для theta
print(rowSums(result2$theta))  # Должно быть вектор единиц длиной D
# Проверка для phi
print(rowSums(result2$phi))    # Должно быть вектор единиц длиной T