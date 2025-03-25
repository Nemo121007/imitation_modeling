# Загрузка необходимых библиотек
library(MASS)  # Для mvrnorm (многомерное нормальное распределение)

# Константы
delta_t <- 0.5
alpha <- 0.6
P0 <- 90
beta <- 3
delta <- 1.5
sigma <- 1  # Предположение среднеквадратичного отклонения
N <- 1000   # Количество частиц

# Определение системных матриц
A_tilde <- matrix(c(1, delta_t, delta_t^2 / 2, 
                    0, 1, delta_t, 
                    0, 0, alpha), nrow=3, byrow=TRUE)
B_tilde <- matrix(c(delta_t^2 / 2, delta_t, 0), nrow=3)
C_tilde <- matrix(c(delta_t^2 / 2, delta_t, 0), nrow=3)
A <- rbind(cbind(A_tilde, matrix(0, 3, 3)), 
           cbind(matrix(0, 3, 3), A_tilde))
B <- rbind(cbind(B_tilde, matrix(0, 3, 1)), 
           cbind(matrix(0, 3, 1), B_tilde))
C <- rbind(cbind(C_tilde, matrix(0, 3, 1)), 
           cbind(matrix(0, 3, 1), C_tilde))

# Состояния управления (уникальные)
u_states <- matrix(c(0, 0, 3.5, 0, 0, 3.5, -3.5, 0, 0, -3.5), 
                   ncol=2, byrow=TRUE)
num_u <- nrow(u_states)

# Матрица переходных вероятностей (пример: сумма по строке = 1)
P <- matrix(c(16, 1, 1, 1, 1,
              1, 16, 1, 1, 1,
              1, 1, 16, 1, 1,
              1, 1, 1, 16, 1,
              1, 1, 1, 1, 16), nrow=5, byrow=TRUE) / 20

# Чтение данных
stations <- read.table("stations.txt", header=FALSE)
stations <- t(stations)
rssi_data <- t(read.table("RSSI-measurements.txt", header=FALSE)) 

L <- ncol(rssi_data)  # Число базовых станций
T <- nrow(rssi_data)  # Число временных шагов

# Инициализация частиц
x_particles <- t(mvrnorm(N, mu=rep(0, 6), 
                         Sigma=diag(c(500, 5, 5, 200, 5, 5))))
u_particles <- u_states[sample(1:num_u, N, replace=TRUE), ]

# Функция для предсказания частиц
predict_particles <- function(x, u) {
  w <- t(mvrnorm(N, mu=c(0, 0), Sigma=sigma^2 * diag(2)))
  A %*% x + B %*% t(u) + C %*% w
}

# Функция для вычисления правдоподобия
compute_likelihood <- function(x, y, stations) {
  d <- matrix(0, nrow=L, ncol=N)
  for (l in 1:L) {
    d[l, ] <- sqrt((x[1, ] - stations[l, 1])^2 + (x[4, ] - stations[l, 2])^2)
  }
  mu <- P0 - 10 * beta * log10(d)
  likelihood <- matrix(0, nrow=L, ncol=N)
  for (l in 1:L) {
    likelihood[l, ] <- dnorm(y[l], mean=mu[l, ], sd=delta, log=TRUE)
  }
  colSums(likelihood)
}

# Функция ресемплинга
systematic_resample <- function(weights) {
  u <- runif(1, 0, 1/N)
  cum_weights <- cumsum(weights)
  indices <- rep(0, N)
  for (i in 1:N) {
    u_i <- u + (i - 1) / N
    indices[i] <- which(cum_weights >= u_i)[1]
  }
  indices
}

# Хранение траектории
x_hat <- matrix(0, nrow=6, ncol=T)

# Основной цикл BPF
for (n in 1:T) {
  # Предсказание
  u_indices <- apply(u_particles, 1, function(u) {
    current_index <- which(apply(u_states, 1, function(state) all(abs(state - u) < 1e-6)))[1]
    sample(1:num_u, 1, prob=P[current_index, ])
  })
  u_particles <- u_states[u_indices, ]
  x_particles <- predict_particles(x_particles, u_particles)
  
  # Обновление весов
  y_n <- as.numeric(rssi_data[n, ])  # Вектор RSSI для текущего временного шага (длина 6)
  log_weights <- compute_likelihood(x_particles, y_n, stations)
  max_log <- max(log_weights)
  weights <- exp(log_weights - max_log)  # Для численной стабильности
  weights <- weights / sum(weights)     # Нормализация
  
  # Ресемплинг
  indices <- systematic_resample(weights)
  x_particles <- x_particles[, indices]
  u_particles <- u_particles[indices, ]
  
  # Оценка состояния
  x_hat[, n] <- rowMeans(x_particles)
}

# Визуализация траектории
plot(x_hat[1, ], x_hat[4, ], type="l", xlab="x", ylab="y", 
     main="Оцененная траектория объекта")
points(stations[, 1], stations[, 2], col="red", pch=16)  # Базовые станции
# Добавляем координаты станций
labels <- paste("(", round(stations[, 1], 1), ",", round(stations[, 2], 1), ")", sep="")
text(stations[, 1], stations[, 2], labels=labels, pos=3, col="red", cex=0.8)

# Определение диапазона осей
x_range <- range(x_hat[1, ], na.rm=TRUE)
y_range <- range(x_hat[4, ], na.rm=TRUE)
x_stations_range <- range(stations[, 1], na.rm=TRUE)
y_stations_range <- range(stations[, 2], na.rm=TRUE)
xlim <- c(min(x_range[1], x_stations_range[1]), max(x_range[2], x_stations_range[2]))
ylim <- c(min(y_range[1], y_stations_range[1]), max(y_range[2], y_stations_range[2]))

plot(x_hat[1, ], x_hat[4, ], type="l", xlab="x", ylab="y", 
     main="Оцененная траектория объекта", xlim=xlim, ylim=ylim)
points(stations[, 1], stations[, 2], col="red", pch=19, cex=1.5)
labels <- paste("(", round(stations[, 1], 1), ",", round(stations[, 2], 1), ")", sep="")
text(stations[, 1], stations[, 2], labels=labels, pos=3, col="red", cex=0.8)
