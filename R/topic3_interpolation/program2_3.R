# Program 2-3: Tabla de diferencias divididas
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 2
# Equivalente a CSL/F2-3.FOR

newton_divided <- function(x_data, f_data, x_target = NULL) {
  n  <- length(x_data)
  dd <- matrix(0, nrow = n, ncol = n)
  dd[, 1] <- f_data

  for (k in 2:n) {
    for (i in 1:(n - k + 1)) {
      dd[i, k] <- (dd[i + 1, k - 1] - dd[i, k - 1]) / (x_data[i + k - 1] - x_data[i])
    }
  }

  # Tabla (mismo formato que el libro)
  cat(" I    X(I)     ", paste(sprintf("F(I,I+%d)  ", 0:(n-1)), collapse=""), "\n")
  for (i in 1:n) {
    j <- n - i
    vals <- dd[i, 1:(j + 1)]
    cat(sprintf("%2d  %7.5f  ", i - 1, x_data[i]))
    cat(sprintf("%9.5f", vals))
    cat("\n")
  }

  if (!is.null(x_target)) {
    result <- dd[1, 1]
    term   <- 1
    for (k in 1:(n - 1)) {
      term   <- term * (x_target - x_data[k])
      result <- result + term * dd[1, k + 1]
    }
    return(result)
  }
}

# ── Datos del libro (cos(x), ejemplo 2-3) ──────────────────────
cat("=== PROGRAMA 2-3: Tabla de Diferencias Divididas ===\n\n")

x <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.2, 1.3)
f <- c(0.99750, 0.99002, 0.96040, 0.88120, 0.76520, 0.67113, 0.62009)

cat("Datos: f(x) = cos(x)\n\n")
res <- newton_divided(x, f, x_target = 0.5)

cat(sprintf("\nInterpolacion en x = 0.5\n"))
cat(sprintf("  Valor interpolado : %.6f\n", res))
cat(sprintf("  Valor exacto      : %.6f\n", cos(0.5)))
cat(sprintf("  Error absoluto    : %.2e\n", abs(cos(0.5) - res)))
