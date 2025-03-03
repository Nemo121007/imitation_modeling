# Устанавливаем начальное значение для воспроизводимости
set.seed(42)

# Функция для генерации одного SAW
generate_SAW <- function(n) {
  path <- list(c(0, 0))
  visited <- list()
  visited[["0,0"]] <- TRUE
  current <- c(0, 0)
  weight <- 1
  
  for (t in 1:n) {
    moves <- list(c(0, 1), c(0, -1), c(1, 0), c(-1, 0))
    valid_moves <- list()
    
    for (move in moves) {
      new_pos <- current + move
      key <- paste(new_pos[1], new_pos[2], sep = ",")
      if (is.null(visited[[key]])) {
        valid_moves[[length(valid_moves) + 1]] <- new_pos
      }
    }
    
    k_t <- length(valid_moves)
    if (k_t == 0) {
      return(NULL)  # Возвращаем NULL при отсутствии ходов
    }
    
    idx <- sample(1:k_t, 1)
    new_pos <- valid_moves[[idx]]
    
    path[[t + 1]] <- new_pos
    visited[[paste(new_pos[1], new_pos[2], sep = ",")]] <- TRUE
    current <- new_pos
    weight <- weight * k_t
  }
  
  final_pos <- path[[n + 1]]
  distance <- sqrt(sum(final_pos^2))
  
  return(list(weight = weight, distance = distance))
}



# Параметры задачи
n <- 20
N <- 10000

# Генерация N корректных путей
results <- list()
while (length(results) < N) {
  result <- generate_SAW(n)
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Извлечение расстояний и весов
distances <- sapply(results, function(x) x$distance)
weights <- sapply(results, function(x) x$weight)

# Оценка среднего расстояния
mu_hat <- sum(distances * weights) / sum(weights)

# Вывод результата
cat("Оценка среднего расстояния:", mu_hat, "\n")