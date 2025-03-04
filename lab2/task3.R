# Параметры для рулетки
Wp <- 0.0001  # порог веса
m <- 10      # коэффициент для рулетки

matrix_size <- 100  # Количество интервалов по каждой оси

# Определяем параметры слоёв
layers <- data.frame(
  layer = 1:5,
  mu_a = c(32, 23, 40, 23, 46),
  # mu_a = c(1e-12, 1e-12, 1e-12, 1e-12, 1e-12),
  mu_s = c(165, 227, 246, 227, 253),
  g    = c(0.72, 0.72, 0.72, 0.72, 0.72),
  # g    = c(0.972, 0.972, 0.972, 0.972, 0.972),
  d    = c(0.01, 0.02, 0.02, 0.09, 0.06)
)

mu_a_sphere <- 51
mu_s_sphere <- 186
g_sphere <- 0.8
l <- 0.02     # Глубина центра шара, см
r <- 0.01    # Радиус шара, см

# Вычисляем границы слоёв по оси Z
layers$z_start <- c(0, head(cumsum(layers$d), -1))
layers$z_end   <- cumsum(layers$d)
max_depth <- 0.2





# Функция для определения параметров слоя по координате z
get_medium_params <- function(pos, layers, l, r) {
  x <- pos["x"]
  y <- pos["y"]
  z <- pos["z"]
  
  # Проверка, находится ли точка внутри шара
  if (x^2 + y^2 + (z - l)^2 <= r^2) {
    # Параметры неоднородности
    return(data.frame(mu_a = mu_a_sphere, mu_s = mu_s_sphere, g = g_sphere))
  } else {
    # Параметры слоя ткани
    layer_params <- get_layer_params(z, layers)
    if (is.null(layer_params)) {
      return(NULL) # Фотон покинул среду
    } else {
      return(layer_params)
    }
  }
}




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
  if (abs(g) < 1e-6) {
    # Для g, близкого к нулю, фазовая функция становится равномерной
    cos_theta <- 2 * xi - 1
  } else {
    term <- (1 - g^2) / (1 - g + 2 * g * xi)
    cos_theta <- (1 + g^2 - term^2) / (2 * g)
    # Ограничиваем значения, чтобы гарантировать, что cos_theta ∈ [-1, 1]
    # cos_theta <- max(min(cos_theta, 1), -1)
  }
  sin_theta <- sqrt(1 - cos_theta^2)
  phi <- runif(1, 0, 2 * pi)
  
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
  # Нормируем вектор 
  u_new <- u_new / sqrt(sum(u_new^2))
  return(u_new)
}





# Функция симуляции одного фотона
simulate_photon <- function(layers, Wp, m, l, r) {
  pos <- c(x = 0, y = 0, z = 0)
  u <- c(0, 0, 1)
  W <- 1
  absorbed_energy <- list()
  
  while(TRUE) {
    current_params <- get_medium_params(pos, layers, l, r)
    if(is.null(current_params)) {
      if(pos["z"] < 0) {
        return(list(status = "reflected", absorption = absorbed_energy))
      } else if(pos["z"] >= max_depth) {
        return(list(status = "transmitted", absorption = absorbed_energy))
      }
    }
    
    mu_a <- current_params$mu_a
    mu_s <- current_params$mu_s
    g    <- current_params$g
    
    # Убедимся, что mu_a и mu_s не слишком малы (меньше машинного эпсилона)
    mu_a <- max(mu_a, 1e-10)  # Минимальное значение для mu_a
    mu_s <- max(mu_s, 1e-10)  # Минимальное значение для mu_s
    
    mu   <- mu_a + mu_s
    
    # Вероятность поглощения
    a1 = runif(1)
    a2 = mu_a / mu
    if (a1 <=  mu_a / mu) {
      absorbed_energy[[length(absorbed_energy) + 1]] <- list(x = pos["x"], z = pos["z"], dW = W)
      return(list(status = "absorbed", absorption = absorbed_energy))
    }
    
    s <- rexp(1, rate = mu)
    pos_new <- pos + s * u
    pos <- pos_new
    
    if(pos["z"] < 0) {
      return(list(status = "reflected", absorption = absorbed_energy))
    } else if(pos["z"] >= max_depth) {
      return(list(status = "transmitted", absorption = absorbed_energy))
    }
    
    dW <- W * (mu_a / mu)
    W <- W - dW
    absorbed_energy[[length(absorbed_energy) + 1]] <- list(x = pos["x"], z = pos["z"], dW = dW)
    
    if(W < Wp) {
      if(runif(1) <= 1/m) {
        W <- m * W
      } else {
        absorbed_energy[[length(absorbed_energy) + 1]] <- list(x = pos["x"], z = pos["z"], dW = W)
        return(list(status = "absorbed", absorption = absorbed_energy))
      }
    }
    
    scatter <- sample_scattering(g)
    u <- update_direction(u, scatter)
  }
}




# Основной цикл симуляции
# set.seed(42)
n_photons <- 10000
n_absorbed <- 0
n_reflected <- 0
n_transmitted <- 0
absorption_all <- list()

for (i in 1:n_photons) {
  res <- simulate_photon(layers, Wp, m, l, r)
  if(res$status == "absorbed") {
    n_absorbed <- n_absorbed + 1
    absorption_all <- c(absorption_all, res$absorption)
  } else if(res$status == "reflected") {
    n_reflected <- n_reflected + 1
  } else if(res$status == "transmitted") {
    n_transmitted <- n_transmitted + 1
  }
}

# Вывод долей
cat("Доля поглощённых фотонов:", n_absorbed / n_photons, "\n")
cat("Доля отражённых фотонов:", n_reflected / n_photons, "\n")
cat("Доля прошедших фотонов:", n_transmitted / n_photons, "\n")

# Построение карты распределения поглощённой энергии
if(length(absorption_all) > 0) {
  abs_df <- do.call(rbind, lapply(absorption_all, function(ev) {
    data.frame(x = ev$x, z = ev$z, dW = ev$dW)
  }))
  
  # Определяем диапазоны координат
  x_range <- range(abs_df$x)
  # Ось Z фиксируется от 0 до max_depth равномерно
  z_range <- c(0, max_depth)
  
  x_bins <- seq(x_range[1], x_range[2], length.out = matrix_size + 1)
  z_bins <- seq(z_range[1], z_range[2], length.out = matrix_size + 1)
  
  energy_matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  
  for(j in 1:nrow(abs_df)) {
    ix <- findInterval(abs_df$x[j], x_bins, rightmost.closed = TRUE)
    iz <- findInterval(abs_df$z[j], z_bins, rightmost.closed = TRUE)
    if(ix >= 1 && ix <= matrix_size && iz >= 1 && iz <= matrix_size) {
      energy_matrix[iz, ix] <- energy_matrix[iz, ix] + abs_df$dW[j]
    }
  }
  
  energy_matrix <- t(energy_matrix)
  energy_matrix <- energy_matrix[, ncol(energy_matrix):1]
  
  eps <- 1e-12
  log_energy_matrix <- log(energy_matrix + eps)
  
  image(x = x_bins, y = z_bins, z = log_energy_matrix,
        xlab = "x (см)", ylab = "z (см)",
        main = "Карта распределения поглощённой энергии (log)",
        col = heat.colors(12))
}