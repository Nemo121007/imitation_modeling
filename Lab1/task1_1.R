set.seed(42)  # Фиксируем сид

# Определяем функцию плотности
# f_x определена на интервале [0, 2π]. В остальном она равна 0
f_x <- function(x) {
  (1 / (4 * pi)) * (2 + cos(x)) * (x >= 0 & x <= 2 * pi)
}

# Генерация выборки методом принятия-отбора
n <- 10000  # Размер выборки
samples <- numeric(n)
count <- 0

while (count < n) {
  x_candidate <- runif(1, 0, 2 * pi)  # Генерируем x из U(0,2π)
  
  # max f_x =  3 / (4 * pi)
  y_candidate <- runif(1, 0, 3 / (4 * pi))  # Верхняя граница
  
  if (y_candidate <= f_x(x_candidate)) {
    count <- count + 1
    samples[count] <- x_candidate
  }
}

# Визуализация полученной выборки
hist(samples, probability = TRUE, breaks = 30, col = "lightblue", main = "Гистограмма выборки и теоретическая плотность")

# Добавляем плотность f_x(x) к гистограмме
curve(f_x(x), add = TRUE, col = "red", lwd = 2)

# Нормирование функции для теста Колмогорова-Смирнова
norm_const <- integrate(f_x, lower = 0, upper = 2 * pi)$value
F_x <- Vectorize(function(x) {
  integrate(f_x, lower = 0, upper = x)$value / norm_const
})

# Проверка гипотезы согласия
ks_test_result <- ks.test(samples, F_x)
print(ks_test_result)