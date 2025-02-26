# Параметры для рулетки
Wp <- 0.001  # порог веса
m <- 10      # коэффициент для рулетки

# Определяем параметры слоёв в виде data.frame
layers <- data.frame(
  layer = 1:5,
  mu_a = c(32, 23, 40, 23, 46),
  mu_s = c(165, 227, 246, 227, 253),
  g    = c(0.72, 0.72, 0.72, 0.72, 0.72),
  d    = c(0.01, 0.02, 0.02, 0.09, 0.06)
)

# Вычисляем границы слоёв по оси Z (кумулятивная сумма)
layers$z_start <- c(0, head(cumsum(layers$d), -1))
layers$z_end   <- cumsum(layers$d)
max_depth <- max(layers$z_end)

# Функция для определения параметров слоя по координате z
get_layer_params <- function(z, layers) {
  idx <- which(z >= layers$z_start & z < layers$z_end)
  if(length(idx) == 0) {
    return(NULL)
  } else {
    return(layers[idx, ])
  }
}

# Функция для генерации угла рассеяния по Хени-Гринштейну
sample_scattering <- function(g) {
  xi <- runif(1)
  if(abs(g) < 1e-6) {
    cos_theta <- 2*xi - 1
  } else {
    temp <- (1 - g^2) / (1 - g + 2*g*xi)
    cos_theta <- (1 + g^2 - temp^2) / (2*g)
    # Ограничим случайные вычисления, если возникают погрешности
    cos_theta <- max(min(cos_theta, 1), -1)
  }
  sin_theta <- sqrt(1 - cos_theta^2)
  phi <- runif(1, 0, 2*pi)
  return(list(cos_theta = cos_theta, sin_theta = sin_theta, phi = phi))
}

# Функция для обновления направления после рассеяния
update_direction <- function(u, scatter) {
  # u - текущий вектор направления (ux, uy, uz)
  cos_theta <- scatter$cos_theta
  sin_theta <- scatter$sin_theta
  phi <- scatter$phi
  
  ux <- u[1]; uy <- u[2]; uz <- u[3]
  
  if(abs(uz) > 0.9999) {
    # Если направление почти вертикальное
    u_new <- c(sin_theta*cos(phi),
               sin_theta*sin(phi),
               sign(uz)*cos_theta)
  } else {
    factor <- sqrt(1 - uz^2)
    u_new <- c(sin_theta*(ux*uz*cos(phi) - uy*sin(phi))/factor + ux*cos_theta,
               sin_theta*(uy*uz*cos(phi) + ux*sin(phi))/factor + uy*cos_theta,
               -sin_theta*cos(phi)*factor + uz*cos_theta)
  }
  # Нормируем вектор (на всякий случай)
  u_new <- u_new / sqrt(sum(u_new^2))
  return(u_new)
}

# Функция симуляции одного фотона
simulate_photon <- function(layers, Wp, m) {
  # Начальные условия
  pos <- c(x = 0, y = 0, z = 0)
  u <- c(0, 0, 1)   # направление вдоль оси Z (перпендикулярно поверхности)
  W <- 1           # начальный вес
  absorbed_energy <- list()  # список для хранения точек поглощения и величины ΔW
  
  while(TRUE) {
    # Определяем текущий слой
    current_layer <- get_layer_params(pos["z"], layers)
    if(is.null(current_layer)) {
      # Если фотон вне ткани (z < 0 или z > max_depth) => завершение симуляции
      return(list(status = "transmitted", absorption = absorbed_energy))
    }
    
    # Извлекаем параметры текущего слоя
    mu_a <- current_layer$mu_a
    mu_s <- current_layer$mu_s
    g    <- current_layer$g
    mu   <- mu_a + mu_s
    
    # Генерируем шаг s из экспоненциального распределения
    s <- rexp(1, rate = mu)
    
    # Обновляем позицию
    pos_new <- pos + s * u
    
    # (Для упрощения не разыгрываем пересечение границ: считаем, что параметр среды постоянен на шаге)
    pos <- pos_new
    
    # Если фотон вышел за пределы ткани
    if(pos["z"] < 0 || pos["z"] >= max_depth) {
      return(list(status = "transmitted", absorption = absorbed_energy))
    }
    
    # Поглощение: вычисляем ΔW и обновляем вес
    dW <- W * (mu_a / mu)
    W <- W - dW
    # Записываем поглощённую энергию и координаты точки (можно брать, например, x и z)
    absorbed_energy[[length(absorbed_energy) + 1]] <- list(x = pos["x"], z = pos["z"], dW = dW)
    
    # Если вес становится очень мал, применяем правило рулетки
    if(W < Wp) {
      if(runif(1) < 1/m) {
        W <- m * W  # фотон выживает с увеличением веса
      } else {
        return(list(status = "absorbed", absorption = absorbed_energy))
      }
    }
    
    # Рассеяние: генерируем углы рассеяния и обновляем направление
    scatter <- sample_scattering(g)
    u <- update_direction(u, scatter)
  }
}

# Основной цикл симуляции
set.seed(123)  # для воспроизводимости
n_photons <- 10000
results <- vector("list", n_photons)
absorption_all <- list()  # для накопления всех событий поглощения
n_absorbed <- 0
n_transmitted <- 0

for(i in 1:n_photons) {
  res <- simulate_photon(layers, Wp, m)
  results[[i]] <- res
  # Обрабатываем список событий поглощения
  if(res$status == "absorbed") {
    n_absorbed <- n_absorbed + 1
    absorption_all <- c(absorption_all, res$absorption)
  } else if(res$status == "transmitted") {
    n_transmitted <- n_transmitted + 1
  }
}

cat("Количество поглощённых фотонов:", n_absorbed, "\n")
cat("Количество прошедших фотонов:", n_transmitted, "\n")
cat("Коэффициент пропускания:", n_transmitted / n_photons, "\n")

# Построение карты распределения поглощённой энергии
# Сначала извлекаем координаты поглощения и соответствующую энергию
if(length(absorption_all) > 0) {
  abs_df <- do.call(rbind, lapply(absorption_all, function(ev) {
    data.frame(x = ev$x, z = ev$z, dW = ev$dW)
  }))
  
  # Определяем диапазоны координат
  x_range <- range(abs_df$x)
  z_range <- range(abs_df$z)
  
  # Разбиваем диапазоны на 10 интервалов
  x_bins <- seq(x_range[1], x_range[2], length.out = 11)
  z_bins <- seq(z_range[1], z_range[2], length.out = 11)
  
  # Создаём матрицу для суммарной энергии
  energy_matrix <- matrix(0, nrow = 10, ncol = 10)
  
  # Аккумулируем энергию в ячейках
  for(j in 1:nrow(abs_df)) {
    ix <- findInterval(abs_df$x[j], x_bins, rightmost.closed = TRUE)
    iz <- findInterval(abs_df$z[j], z_bins, rightmost.closed = TRUE)
    # Если значение попало в 11-ю границу, уменьшаем индекс до 10
    if(ix == 11) ix <- 10
    if(iz == 11) iz <- 10
    energy_matrix[iz, ix] <- energy_matrix[iz, ix] + abs_df$dW[j]
  }
  
  energy_matrix <- t(energy_matrix)
  
  eps <- 1e-12
  log_energy_matrix <- log(energy_matrix + eps)
  
  # Визуализация с помощью логарифмической шкалы
  image(x = x_bins, y = z_bins, z = log_energy_matrix,
        xlab = "x (см)", ylab = "z (см)",
        main = "Карта распределения поглощённой энергии (log)",
        col = heat.colors(12))
}

