library(deSolve)

# Чтение данных из файла data.txt (только значения y, без заголовка)
y <- scan("data.txt", what = numeric())  # Считываем как числовой вектор
n <- length(y)

# Генерация вектора времени
t <- seq(0, (n - 1) * 0.001, by = 0.001)  # От 0 до 2.998 с шагом 0.001

# Параметры задачи
CA0 <- 10    # Начальная концентрация A
CB0 <- 15    # Начальная концентрация B
delta <- 0.2 # Стандартное отклонение шума

# Модель дифференциальных уравнений
reaction_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dCA <- -k1 * CA + k2 * CB
    dCB <- k1 * CA - k2 * CB
    list(c(dCA, dCB))
  })
}

# Функция для решения системы и получения C_A(t) с проверками
solve_CA <- function(k1, k2, t) {
  state <- c(CA = CA0, CB = CB0)
  parameters <- c(k1 = k1, k2 = k2)
  out <- tryCatch(
    ode(y = state, times = t, func = reaction_model, parms = parameters, hmax = 0.001),
    error = function(e) NULL
  )
  if (is.null(out) || any(is.na(out[, "CA"]))) {
    return(rep(NA, length(t)))  # Возвращаем NA, если решение не удалось
  }
  return(out[, "CA"])
}

# Логарифм правдоподобия 
log_likelihood <- function(k1, k2, y, t, delta) {
  CA_pred <- solve_CA(k1, k2, t)
  ll <- sum(dnorm(y, mean = CA_pred, sd = delta, log = TRUE))
  return(ll)
}

# Логарифм априорного распределения
log_prior <- function(k1, k2) {
  if (k1 > 0 && k1 < 10 && k2 > 0 && k2 < 10) {
    return(-2 * log(10)) # log(1/100)
  } else {
    return(-Inf)
  }
}

# Логарифм апостериорного распределения (π*(x) ∝ π(x))
log_posterior <- function(k1, k2, y, t, delta) {
  ll <- log_likelihood(k1, k2, y, t, delta)
  prior <- log_prior(k1, k2)
  post <- ll + prior
  return(post)
}

# Алгоритм Метрополиса-Гастингса
metropolis_hastings <- function(y, t, delta, n_iter, k1_init, k2_init, sigma_prop) {
  k1_chain <- numeric(n_iter)
  k2_chain <- numeric(n_iter)
  
  # 1. Выбираем начальное состояние X₀ = (k1_init, k2_init)
  k1_current <- k1_init
  k2_current <- k2_init
  log_post_current <- log_posterior(k1_current, k2_current, y, t, delta)
  
  # 2. Для n = 0, 1, ..., n_iter - 1 выполняем итерации
  for (i in 1:n_iter) {
    # 3. Генерируем кандидата Y ~ q(y|Xₙ), где q(y|x) = N(x, σ²)
    k1_prop <- rnorm(1, k1_current, sigma_prop)
    k2_prop <- rnorm(1, k2_current, sigma_prop)
    
    # 4. Вычисляем вероятность принятия α(Xₙ, Y)
    # α(Xₙ, Y) = min{ [π(Y) q(Xₙ|Y)] / [π(Xₙ) q(Y|Xₙ)], 1 }, но q симметрично, поэтому α = min{ π(Y)/π(Xₙ), 1 }
    log_post_prop <- log_posterior(k1_prop, k2_prop, y, t, delta)
    alpha <- min(exp(log_post_prop - log_post_current), 1)
    
    # 5. С вероятностью α принимаем Y, иначе остаемся в Xₙ
    if (runif(1) < alpha) {
      k1_current <- k1_prop
      k2_current <- k2_prop
      log_post_current <- log_post_prop
    }
    
    # Сохраняем текущее состояние в цепь
    k1_chain[i] <- k1_current
    k2_chain[i] <- k2_current
  }
  
  return(list(k1_chain = k1_chain, k2_chain = k2_chain))
}

# Параметры алгоритма
n_iter <- 100      # Количество итераций
k1_init <- 5        # Начальное значение k1
k2_init <- 7        # Начальное значение k2
sigma_prop <- 0.05  # Стандартное отклонение σ для q(y|x)

# Запуск алгоритма
set.seed(42) # Для воспроизводимости
result <- metropolis_hastings(y, t, delta, n_iter, k1_init, k2_init, sigma_prop)

# Анализ результатов
burn_in <- 0.1 * n_iter # Период прогрева
k1_chain <- result$k1_chain[(burn_in + 1):n_iter]
k2_chain <- result$k2_chain[(burn_in + 1):n_iter]

# Оценки параметров
k1_est <- mean(k1_chain)
k2_est <- mean(k2_chain)
k1_sd <- sd(k1_chain)
k2_sd <- sd(k2_chain)

# Вывод результатов
cat("Оценка k1:", k1_est, "(стандартное отклонение:", k1_sd, ")\n")
cat("Оценка k2:", k2_est, "(стандартное отклонение:", k2_sd, ")\n")

# Визуализация цепей
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(result$k1_chain, type = "l", main = "Цепь Маркова для k1", ylab = "k1")
abline(h = k1_est, col = "red")
plot(result$k2_chain, type = "l", main = "Цепь Маркова для k2", ylab = "k2")
abline(h = k2_est, col = "red")

# Гистограммы апостериорных распределений
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(k1_chain, breaks = 50, main = "Апостериорное распределение k1", xlab = "k1")
abline(v = k1_est, col = "red")
hist(k2_chain, breaks = 50, main = "Апостериорное распределение k2", xlab = "k2")
abline(v = k2_est, col = "red")

# Вывод результатов
cat("Оценка k1:", k1_est, "(стандартное отклонение:", k1_sd, ")\n")
cat("Оценка k2:", k2_est, "(стандартное отклонение:", k2_sd, ")\n")