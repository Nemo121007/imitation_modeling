set.seed(42)

simulate_SAW_classic <- function(n, N) {
  # n = длина блуждания
  # N = число генерируемых траекторий
  
  # Функция для определения всех доступных ходов (4 направления),
  # исключая уже посещённые точки
  possible_moves <- function(x, visited) {
    moves <- list(
      c(x[1] + 1, x[2]),  # вправо
      c(x[1] - 1, x[2]),  # влево
      c(x[1], x[2] + 1),  # вверх
      c(x[1], x[2] - 1)   # вниз
    )
    valid_moves <- list()
    for (m in moves) {
      if (!paste(m, collapse=",") %in% visited) {
        valid_moves[[length(valid_moves) + 1]] <- m
      }
    }
    return(valid_moves)
  }
  
  distances <- numeric(N)  # финальные расстояния
  
  for (i in seq_len(N)) {
    path <- list(c(0,0))
    visited <- c("0,0")
    stuck <- FALSE
    
    for (step_idx in seq_len(n)) {
      current_pos <- path[[step_idx]]
      valid_moves <- possible_moves(current_pos, visited)
      m_k <- length(valid_moves)
      
      if (m_k == 0) {
        # Заперлось — досрочно завершим блуждание
        stuck <- TRUE
        break
      }
      
      # Равномерно выбираем из доступных ходов
      chosen_idx <- sample.int(m_k, size=1)
      chosen_move <- valid_moves[[chosen_idx]]
      
      path[[step_idx + 1]] <- chosen_move
      visited <- c(visited, paste(chosen_move, collapse=","))
    }
    
    # Финальная позиция (даже если шагов меньше n)
    final_pos <- path[[length(path)]]
    dist_final <- sqrt(final_pos[1]^2 + final_pos[2]^2)
    distances[i] <- dist_final
  }
  
  # Оценка мат. ожидания (среднего расстояния)
  mean_dist <- mean(distances)
  # Оценка дисперсии (классическое var() в R — несмещённая оценка)
  var_dist <- var(distances)
  
  return(list(mean = mean_dist, variance = var_dist))
}

simulate_SAW_importance <- function(n, N, alpha) {
  # n     = длина блуждания
  # N     = число генерируемых траекторий
  # alpha = параметр смещения для важностной выборки
  
  possible_moves <- function(x, visited) {
    moves <- list(
      c(x[1] + 1, x[2]),
      c(x[1] - 1, x[2]),
      c(x[1], x[2] + 1),
      c(x[1], x[2] - 1)
    )
    valid_moves <- list()
    for (m in moves) {
      if (!paste(m, collapse=",") %in% visited) {
        valid_moves[[length(valid_moves) + 1]] <- m
      }
    }
    return(valid_moves)
  }
  
  # Накопление:
  # - w_i * X_i (для расчёта E[X])
  # - w_i * X_i^2 (для расчёта E[X^2])
  # - w_i (для нормировки)
  sum_w_x   <- 0
  sum_w_x2  <- 0
  sum_w     <- 0
  
  for (i in seq_len(N)) {
    path <- list(c(0,0))
    visited <- c("0,0")
    
    # Логарифмы в целях меньшего изменения порядков значени
    log_f <- 0  # log(f(x_i))
    log_g <- 0  # log(g(x_i))
    
    for (step_idx in seq_len(n)) {
      current_pos <- path[[step_idx]]
      valid_moves <- possible_moves(current_pos, visited)
      m_k <- length(valid_moves)
      
      if (m_k == 0) {
        # Заперлось
        break
      }
      
      # Вероятность при f: равномерно среди m_k => p_f = 1/m_k
      # Вероятность при g: пропорционально exp(alpha * distance)
      g_weights <- numeric(m_k)
      for (j in seq_len(m_k)) {
        new_pos <- valid_moves[[j]]
        dist_new <- sqrt(new_pos[1]^2 + new_pos[2]^2)
        g_weights[j] <- exp(alpha * dist_new)
      }
      g_weights <- g_weights / sum(g_weights)
      
      # Выбираем шаг из valid_moves по g_weights
      chosen_idx <- sample.int(m_k, size=1, prob=g_weights)
      
      # Обновляем логи
      # log_f += log(1/m_k)
      log_f <- log_f + (-log(m_k))
      # log_g += log(g_weights[chosen_idx])
      log_g <- log_g + log(g_weights[chosen_idx])
      
      chosen_move <- valid_moves[[chosen_idx]]
      path[[step_idx + 1]] <- chosen_move
      visited <- c(visited, paste(chosen_move, collapse=","))
    }
    
    final_pos <- path[[length(path)]]
    dist_final <- sqrt(final_pos[1]^2 + final_pos[2]^2)
    
    # Вес w_i = f(x_i) / g(x_i) = exp(log_f - log_g)
    w_i <- exp(log_f - log_g)
    
    # Накопим суммы
    sum_w_x  <- sum_w_x  + w_i * dist_final
    sum_w_x2 <- sum_w_x2 + w_i * dist_final^2
    sum_w    <- sum_w    + w_i
  }
  
  # Оценка E[X] = (sum(w_i * X_i)) / (sum w_i)
  E_x <- sum_w_x / sum_w
  # Оценка E[X^2] = (sum(w_i * X_i^2)) / (sum w_i)
  E_x2 <- sum_w_x2 / sum_w
  # Дисперсия = E[X^2] - (E[X])^2
  var_x <- E_x2 - E_x^2
  
  return(list(mean = E_x, variance = var_x))
}



n <- 20     # длина блуждания
N <- 10000  # число траекторий
alpha <- 5

# 1) Классический метод Монте-Карло
res_classic <- simulate_SAW_classic(n, N)

# 2) Важностная выборка (IS)
res_importance <- simulate_SAW_importance(n, N, alpha)

cat("==== РЕЗУЛЬТАТЫ ====\n")
cat("Классический метод Монте-Карло:\n")
cat("  Оценка среднего расстояния =", res_classic$mean, "\n")
cat("  Оценка дисперсии           =", res_classic$variance, "\n\n")

cat("Метод важностной выборки (alpha =", alpha, "):\n")
cat("  Оценка среднего расстояния =", res_importance$mean, "\n")
cat("  Оценка дисперсии           =", res_importance$variance, "\n")
