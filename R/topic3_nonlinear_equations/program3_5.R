# Program 3-5: Método de Newton-Raphson
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 3
# Equivalente a CSL/F3-5.FOR

# ── Función y su derivada (ejemplo 3.5 del libro) ─────────────────────────
# f(x)  = x³ - 5x² + 6x      raíces en 0, 2, 3
# f'(x) = 3x² - 10x + 6
func_newton <- function(x) {
  list(
    y  = x^3 - 5 * x^2 + 6 * x,
    yd = 3 * x^2 - 10 * x + 6
  )
}

# ── Newton-Raphson ─────────────────────────────────────────────────────────
# f_and_deriv: función que devuelve list(y, yd)
# x0:          estimación inicial
# tol:         tolerancia de convergencia  |x_new - x_old| < tol
newton <- function(f_and_deriv, x0, tol) {
  cat(sprintf("\n ESTIMACION INICIAL: %.6f\n", x0))
  cat(sprintf(" %-8s  %-14s  %-14s  %-14s\n", "IT.NO. N", "X(N-1)", "Y(N-1)", "X(N)"))

  x  <- x0
  xb <- x0
  i  <- 0L

  repeat {
    i   <- i + 1L
    res <- f_and_deriv(x)
    x_new <- x - res$y / res$yd

    cat(sprintf(" %5d     %14.6e  %14.6e  %14.6e\n", i, xb, res$y, x_new))

    if (abs(x_new - xb) < tol) {
      cat("-----------------------------------------\n")
      cat(sprintf("   SOLUCION FINAL = %.6f\n", x_new))
      cat("-----------------------------------------\n")
      return(invisible(x_new))
    }

    xb <- x_new
    x  <- x_new
  }
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 3-5: ESQUEMA DE NEWTON ===\n")
cat("    f(x) = x^3 - 5x^2 + 6x     (raices en 0, 2, 3)\n")
cat("    f'(x) = 3x^2 - 10x + 6\n")

tol <- 0.00001

# Caso 1: x0 = 4.0  -> raíz en 3
newton(func_newton, x0 = 4.0, tol)

# Caso 2: x0 = 1.4  -> raíz en 2
newton(func_newton, x0 = 1.4, tol)

# Caso 4: x0 = 0.1  -> raíz en 0
newton(func_newton, x0 = 0.1, tol)
