# Устанавливаем зерно для воспроизводимости
set.seed(42)

# Задаём параметры
n <- 5          # Число злокачественных клеток
m <- 10        # Число здоровых клеток
k <- 3          # Порог для N
N_sim <- 10000  # Число симуляций

# Параметры экспоненциальных распределений
lambda_mal <- rep(2, n)        # λ для злокачественных клеток
lambda_healthy <- rep(1, m)  # λ для здоровых клеток



# **Прямой метод Монте-Карло**
indicator_direct <- numeric(N_sim)
for (i in 1:N_sim) {
  # Генерируем времена гибели злокачественных клеток
  T_mal <- rexp(n, rate = lambda_mal)
  M <- max(T_mal)
  # Генерируем времена гибели здоровых клеток
  T_healthy <- rexp(m, rate = lambda_healthy)
  # Находим k-ое по убыванию время
  T_k <- sort(T_healthy, decreasing = TRUE)[k]
  # Индикатор условия M < T^{(k)}
  indicator_direct[i] <- ifelse(M < T_k, 1, 0)
}
# Оценка вероятности
p_hat_direct <- mean(indicator_direct)
# Дисперсия оценки
var_direct <- var(indicator_direct) / N_sim



# **Условный метод Монте-Карло**
# Генерируем времена гибели здоровых клеток для всех симуляций
T_healthy_matrix <- matrix(rexp(N_sim * m, rate = lambda_healthy), nrow = N_sim)
# Находим T^{(k)} для каждой симуляции
T_k_vec <- apply(T_healthy_matrix, 1, function(x) sort(x, decreasing = TRUE)[k])
# Вычисляем Z = P(M < T^{(k)} | T^{(k)})
lambda_M <- sum(lambda_mal)
Z <- sapply(T_k_vec, function(t) prod(1 - exp(-lambda_mal * t)))
# Оценка вероятности
p_hat_conditional <- mean(Z)
# Дисперсия оценки
var_conditional <- var(Z) / N_sim


# Вывод результатов
list(
  "Прямой метод Монте-Карло:",
  "Оценка P(N >= k):" = p_hat_direct,
  "Дисперсия оценки:" = var_direct,
  "Условный метод Монте-Карло",
  "Оценка P(N >= k):" = p_hat_conditional,
  "Дисперсия оценки:" = var_conditional,
  "Разница в оценках P(N >= k):" = p_hat_direct - p_hat_conditional,
  "Отношение дисперсий (прямой / условный):" = var_direct / var_conditional
)
