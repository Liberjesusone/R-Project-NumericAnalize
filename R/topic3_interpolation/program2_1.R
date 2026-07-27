# Program 2-1: Interpolacion de Lagrange
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 2
# Equivalente a CSL/F2-1.FOR

lagrange <- function(x_data, f_data, xa) {
  n    <- length(x_data)
  yres <- 0

  if (xa < min(x_data) || xa > max(x_data))
    cat("  ADVERTENCIA: X esta en el rango de extrapolacion\n")

  for (i in 1:n) {
    z <- 1.0
    for (j in 1:n) {
      if (i != j) z <- z * (xa - x_data[j]) / (x_data[i] - x_data[j])
    }
    yres <- yres + z * f_data[i]
  }
  return(yres)
}

# ── Datos del libro (ejemplo 2-1) ──────────────────────────────
cat("=== PROGRAMA 2-1: Interpolacion de Lagrange ===\n\n")

x <- c(1,    2,    3,    4   )
f <- c(0.671, 0.620, 0.567, 0.512)

cat("Tabla de valores:\n")
cat(sprintf("  %3s  %8s  %8s\n", "I", "X(I)", "F(I)"))
for (i in seq_along(x))
  cat(sprintf("  %3d  %8.3f  %8.5f\n", i - 1, x[i], f[i]))

# Mismos valores de prueba que el ejemplo de salida del libro
puntos <- c(3.66, 4.5, 0.1)
cat("\n")
for (xa in puntos) {
  cat(sprintf("DAR X = %.2f\n", xa))
  res <- lagrange(x, f, xa)
  cat(sprintf("  RESULTADO: G(%12.5e) = %12.5e\n\n", xa, res))
}
