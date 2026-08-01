# Spline Cúbico Natural
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 2
# Implementacion desde cero — sin usar splinefun() de R

# ── Algoritmo de Thomas: resuelve sistema tridiagonal en O(n) ──────────────
# Es el método más eficiente para este tipo de sistema.
# lower: sub-diagonal   (n-1 elementos)
# main:  diagonal princ (n   elementos)
# upper: super-diagonal (n-1 elementos)
thomas_solve <- function(lower, main, upper, rhs) {
  n <- length(rhs)
  if (n == 0) return(numeric(0))
  b <- main
  d <- rhs
  if (n >= 2) {
    for (i in 2:n) {
      m    <- lower[i - 1] / b[i - 1]
      b[i] <- b[i] - m * upper[i - 1]
      d[i] <- d[i] - m * d[i - 1]
    }
  }
  x <- numeric(n)
  x[n] <- d[n] / b[n]
  if (n >= 2) {
    for (i in (n - 1):1)
      x[i] <- (d[i] - upper[i] * x[i + 1]) / b[i]
  }
  return(x)
}

# ── Ajuste del spline cúbico natural ───────────────────────────────────────
# Devuelve los coeficientes a, b, c, d de cada tramo y los M_i = S''(x_i).
# S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
spline_fit <- function(x_data, f_data) {
  n <- length(x_data) - 1   # número de intervalos
  h <- diff(x_data)          # h[k] = x_{k+1} - x_k

  # Lado derecho del sistema para los M_i interiores
  rhs <- numeric(n - 1)
  for (k in seq_len(n - 1)) {
    rhs[k] <- 6 * ((f_data[k + 2] - f_data[k + 1]) / h[k + 1] -
                    (f_data[k + 1] - f_data[k])   / h[k])
  }

  # Diagonales del sistema tridiagonal
  d_main <- 2 * (h[seq_len(n - 1)] + h[seq_len(n - 1) + 1])
  d_off  <- if (n > 2) h[2:(n - 1)] else numeric(0)

  # Resolver: M_0 = M_n = 0 (condicion de frontera natural)
  M_int <- thomas_solve(d_off, d_main, d_off, rhs)
  M <- c(0, M_int, 0)

  # Coeficientes de cada tramo
  idx <- seq_len(n)
  a   <- f_data[idx]
  b   <- (f_data[idx + 1] - f_data[idx]) / h[idx] - h[idx] * (2*M[idx] + M[idx + 1]) / 6
  cc  <- M[idx] / 2
  d   <- (M[idx + 1] - M[idx]) / (6 * h[idx])

  list(x=x_data, a=a, b=b, c=cc, d=d, M=M, n=n)
}

# ── Evaluación del spline en un vector de puntos ───────────────────────────
spline_eval <- function(sp, x_eval) {
  i  <- findInterval(x_eval, sp$x, rightmost.closed = TRUE)
  i  <- pmax(1L, pmin(i, sp$n))
  dx <- x_eval - sp$x[i]
  sp$a[i] + sp$b[i] * dx + sp$c[i] * dx^2 + sp$d[i] * dx^3
}

# ══════════════════════════════════════════════════════════════════════════════
# DEMO: f(x) = x·sin(x) en [0, π/2] — los mismos 6 nodos del quiz
# ══════════════════════════════════════════════════════════════════════════════
cat("=== SPLINE CÚBICO NATURAL: f(x) = x·sin(x) ===\n\n")

x <- seq(0, pi / 2, length.out = 6)   # h = π/10
f <- x * sin(x)

cat("Datos de entrada:\n")
cat(sprintf("  %2s  %10s  %10s\n", "i", "x_i", "f_i"))
for (i in seq_along(x))
  cat(sprintf("  %2d  %10.6f  %10.6f\n", i - 1L, x[i], f[i]))

# ── Ajustar ───────────────────────────────────────────────────────────────
sp <- spline_fit(x, f)

cat("\nSegundas derivadas en los nodos M_i = S''(x_i):\n")
for (i in seq_along(sp$M))
  cat(sprintf("  M_%d = %10.6f\n", i - 1L, sp$M[i]))

cat("\nCoeficientes de cada tramo  S_i(x) = a + b·dx + c·dx² + d·dx³:\n")
cat(sprintf("  %2s  %11s  %11s  %11s  %11s\n", "i", "a", "b", "c", "d"))
for (i in seq_len(sp$n))
  cat(sprintf("  %2d  %11.6f  %11.6f  %11.6f  %11.6f\n",
              i - 1L, sp$a[i], sp$b[i], sp$c[i], sp$d[i]))

# ── Interpolar en π/4 ─────────────────────────────────────────────────────
xa     <- pi / 4
exacto <- xa * sin(xa)
val_sp <- spline_eval(sp, xa)

cat(sprintf("\nInterpolación en x = π/4 = %.6f\n", xa))
cat(sprintf("  Valor spline   : %.6f\n", val_sp))
cat(sprintf("  Valor exacto   : %.6f\n", exacto))
cat(sprintf("  Error absoluto : %.2e\n", abs(exacto - val_sp)))
cat(sprintf("  Error relativo : %.5f%%\n", abs(exacto - val_sp) / exacto * 100))

# ── Comparación con los métodos del quiz ──────────────────────────────────
cat("\nComparación para x = π/4:\n")
cat(sprintf("  %-22s %10s  %10s\n", "Método", "P(π/4)", "Error rel."))
cat(sprintf("  %-22s %10s  %10s\n", "------", "------", "----------"))
cat(sprintf("  %-22s %10.6f  %10s\n", "2 nodos (lineal)",   0.242703, "56.3%"))
cat(sprintf("  %-22s %10.6f  %10s\n", "4 nodos (Newton)",   0.554170, "0.21%"))
cat(sprintf("  %-22s %10.6f  %9.5f%%\n", "Spline cúbico (6 pt)", val_sp,
            abs(exacto - val_sp)/exacto * 100))
cat(sprintf("  %-22s %10.6f  %10s\n", "Exacto",             exacto,  "—"))

# ── Calidad en todo el intervalo ──────────────────────────────────────────
cat("\nCalidad del spline en todo [0, π/2]:\n")
cat(sprintf("  %10s  %10s  %10s  %10s\n", "x", "Spline", "Exacto", "Error"))
x_test <- seq(0, pi / 2, by = pi / 20)
for (xi in x_test) {
  fi_e <- xi * sin(xi)
  fi_s <- spline_eval(sp, xi)
  cat(sprintf("  %10.6f  %10.6f  %10.6f  %9.2e\n", xi, fi_s, fi_e, abs(fi_e - fi_s)))
}

# ── Verificación contra la función nativa de R ────────────────────────────
cat("\nVerificación contra splinefun() nativo de R:\n")
sp_r <- splinefun(x, f, method = "natural")
cat(sprintf("  Nuestra impl.  : %.8f\n", val_sp))
cat(sprintf("  splinefun()    : %.8f\n", sp_r(xa)))
cat(sprintf("  Diferencia     : %.2e\n", abs(val_sp - sp_r(xa))))
