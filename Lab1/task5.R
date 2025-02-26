library(ggplot2)
library(plotly)
library(dplyr)

# Количество точек
n_points <- 5000
u_max <- 2 * pi  # Задаем максимальный радиус
M <- sqrt(u_max^2 + 1)  # Верхняя граница функции отбора

# Генерация сразу всех значений v
v <- runif(n_points, 0, 2 * pi)

# Генерация u методом Неймана
accepted_u <- numeric(n_points)
i <- 1

while (i <= n_points) {
  u_prime <- runif(1, 0, u_max)  # Кандидат на u
  U <- runif(1, 0, 1)  # Случайная точка для принятия/отбора
  
  if (U <= sqrt(u_prime^2 + 1) / M) {
    accepted_u[i] <- u_prime
    i <- i + 1
  }
}

# Функции параметризации геликоида
x <- accepted_u * cos(v)
y <- accepted_u * sin(v)
z <- v

# Создание датафрейма
data <- data.frame(x = x, y = y, z = z)

# 3D-график точек
p1 <- plot_ly(data = data, x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "markers",
              marker = list(size = 2, color = data$z, colorscale = "Viridis")) %>%
  layout(title = "Равномерные точки на геликоиде")

# Вывод графика
p1


# Генерация сравниваемого распределения
u_test <- sqrt(runif(n_points)) * 2 * pi

# Проведение теста Колмогорова-Смирнова
ks_test_result <- ks.test(accepted_u, u_test)

# Вывод p-значения теста Колмогорова-Смирнова
print(paste("p-value:", ks_test_result))
