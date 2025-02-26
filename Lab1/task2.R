library(ggplot2)
set.seed(42)



# Функция для проверки, находится ли точка внутри треугольника
is_point_in_triangle <- function(pt, v1, v2, v3) {
  # Вспомогательная функция для вычисления знака площади параллелограмма
  sign <- function(p1, p2, p3) {
    return ((p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2]))
  }
  
  d1 <- sign(pt, v1, v2)
  d2 <- sign(pt, v2, v3)
  d3 <- sign(pt, v3, v1)
  
  # Проверяем, имеют ли все три знака одинаковое направление
  has_neg <- (d1 <= 0 && d2 <= 0 && d3 <= 0)
  has_pos <- (d1 >= 0 && d2 >= 0 && d3 >= 0)
  
  return (has_neg || has_pos)
}





# Функция для визуализации треугольников и точек
plot_intersection <- function(tri1, tri2, points, intersect_points) {
  tri1_df <- data.frame(x = c(tri1[,1], tri1[1,1]), y = c(tri1[,2], tri1[1,2]))
  tri2_df <- data.frame(x = c(tri2[,1], tri2[1,1]), y = c(tri2[,2], tri2[1,2]))
  points_df <- data.frame(x = points$x, y = points$y)
  intersect_df <- data.frame(x = intersect_points$x, y = intersect_points$y)
  
  p <- ggplot() +
    # Рисуем первый треугольник
    geom_polygon(data = tri1_df, aes(x = x, y = y), fill = "blue", alpha = 0.4) +
    # Рисуем второй треугольник
    geom_polygon(data = tri2_df, aes(x = x, y = y), fill = "red", alpha = 0.4) +
    # Рисуем случайные точки
    geom_point(data = points_df, aes(x = x, y = y), color = "black", alpha = 0.2) +
    # Рисуем точки, которые попали в пересечение
    geom_point(data = intersect_df, aes(x = x, y = y), color = "green", alpha = 0.6) +
    # Настроим оси и график
    labs(title = "Пересечение треугольников и случайные точки",
         x = "X", y = "Y") +
    theme_minimal()
  
  # Используем print для явного вывода графика
  print(p)
}





# Функция для вычисления пересечения и визуализации
monte_carlo_intersection_with_plot <- function(tri1, tri2, n_points = 10000) {
  # Минимальные и максимальные значения для прямоугольника, содержащего оба треугольника
  min_x <- min(c(tri1[,1], tri2[,1]))
  max_x <- max(c(tri1[,1], tri2[,1]))
  min_y <- min(c(tri1[,2], tri2[,2]))
  max_y <- max(c(tri1[,2], tri2[,2]))
  
  # Генерация случайных точек внутри этого прямоугольника
  points <- data.frame(x = runif(n_points, min_x, max_x), 
                       y = runif(n_points, min_y, max_y))
  
  # Найдем точки, попавшие в оба треугольника (пересечение)
  intersect_points <- points[apply(points, 1, function(pt) 
    is_point_in_triangle(c(pt[1], pt[2]), tri1[1,], tri1[2,], tri1[3,]) && 
      is_point_in_triangle(c(pt[1], pt[2]), tri2[1,], tri2[2,], tri2[3,])), ]
  
  # Площадь пересечения
  count_intersection <- nrow(intersect_points)
  
  # Площадь прямоугольника
  area_rectangle <- (max_x - min_x) * (max_y - min_y)
  
  # Площадь пересечения
  intersection_area <- (count_intersection / n_points) * area_rectangle
  
  # Визуализация
  plot_intersection(tri1, tri2, points, intersect_points)
  
  return(intersection_area)
}





# Пример: Большой треугольник и вложенный в него маленький
tri1 <- matrix(c(0, 0, 5, 0, 0, 5), ncol = 2, byrow = TRUE)
tri2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

# Пример: Два пересекающихся треугольника
tri1 <- matrix(c(1, 1, 5, 1, 3, 4), ncol = 2, byrow = TRUE)  
tri2 <- matrix(c(2, 0, 6, 0, 4, 3), ncol = 2, byrow = TRUE)  

# Пример: Треугольники с одной общей вершиной
# tri1 <- matrix(c(0, 0, 4, 0, 2, 3), ncol = 2, byrow = TRUE)  
# tri2 <- matrix(c(0, 0, 3, -2, 1, -4), ncol = 2, byrow = TRUE)

# Пример: Треугольник, находящийся внутри другого, но не касающийся вершин
# tri1 <- matrix(c(-3, -3, -1, -3, -2, -1), ncol = 2, byrow = TRUE)
# tri2 <- matrix(c(2, 2, 4, 2, 3, 4), ncol = 2, byrow = TRUE)

# Вызов функции Монте-Карло с визуализацией
area_intersection <- monte_carlo_intersection_with_plot(tri1, tri2)
print(paste("Площадь пересечения: ", area_intersection))
