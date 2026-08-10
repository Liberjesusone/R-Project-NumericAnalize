# Program 4-2: Fórmulas cerradas de Newton-Cotes
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 4
# Equivalente a CSL/F4-2.FOR

# ── Tabla 4.2: pesos w_i y alfa = Q/R para cada orden N ───────────────────
# La fórmula es:  I = alfa * h * sum( w_i * f(x_i) ),   h = (b-a)/N
# Cada entrada tiene N+1 pesos. Q y R se guardan aparte.
NC_W <- list(
  "1"  = c(1, 1),
  "2"  = c(1, 4, 1),
  "3"  = c(1, 3, 3, 1),
  "4"  = c(7, 32, 12, 32, 7),
  "5"  = c(19, 75, 50, 50, 75, 19),
  "6"  = c(41, 216, 27, 272, 27, 216, 41),
  "7"  = c(751, 3577, 1323, 2989, 2989, 1323, 3577, 751),
  "8"  = c(989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989),
  "9"  = c(2857, 15741, 1080, 19344, 5778, 5778, 19344, 1080, 15741, 2857),
  "10" = c(16067, 106300, -48525, 272400, -260550, 427368,
           -260550, 272400, -48525, 106300, 16067)
)

NC_Q <- c("1"=1, "2"=1, "3"=3, "4"=2, "5"=5,
          "6"=1, "7"=7, "8"=4, "9"=9, "10"=5)
NC_R <- c("1"=2, "2"=3, "3"=8, "4"=45, "5"=288,
          "6"=140, "7"=17280, "8"=14175, "9"=89600, "10"=299376)

# ── Verificación de consistencia de la tabla ──────────────────────────────
# Una fórmula de Newton-Cotes debe integrar exactamente f(x) = 1.
# Con f = 1:  I = alfa*h*sum(w) = alfa*((b-a)/N)*sum(w) = b-a
# Por lo tanto se debe cumplir:  alfa * sum(w) = N
check_weights <- function() {
  cat("--- Verificacion de la tabla 4.2:  alfa * sum(w) = N ? ---\n")
  cat(sprintf(" %4s  %12s  %14s  %10s  %s\n", "N", "sum(w)", "alfa", "alfa*sum", "OK"))
  for (n in 1:10) {
    k  <- as.character(n)
    sw <- sum(NC_W[[k]])
    al <- NC_Q[k] / NC_R[k]
    pr <- al * sw
    cat(sprintf(" %4d  %12.0f  %14.8f  %10.4f  %s\n",
                n, sw, al, pr, ifelse(abs(pr - n) < 1e-9, "SI", "NO <-- ERROR")))
  }
  cat("\n")
}

# ── Fórmula cerrada de Newton-Cotes de orden N ────────────────────────────
# k = número de datos (puntos de la retícula);  N = k - 1 = orden
newton_cotes <- function(f, a, b, k, verbose = TRUE) {
  n  <- k - 1
  ky <- as.character(n)
  if (is.null(NC_W[[ky]])) stop("Orden fuera de rango (2 <= k <= 11)")

  w  <- NC_W[[ky]]
  q  <- NC_Q[ky]; r <- NC_R[ky]
  al <- q / r                 # alfa
  h  <- (b - a) / n           # espaciamiento de la retícula

  if (verbose) {
    cat(sprintf("\n Q = %14.6e    R = %14.6e    alfa = Q/R = %.8f\n", q, r, al))
    cat("--------------------------------------------------\n")
    cat(sprintf(" %5s  %14s  %14s  %12s\n", "N", "X", "F(X)", "W"))
    cat("--------------------------------------------------\n")
  }

  ans <- 0
  for (j in 0:n) {
    x  <- a + j * h
    fx <- f(x)
    if (verbose) cat(sprintf(" %5d  %14.6e  %14.6e  %12.4f\n", j, x, fx, w[j + 1]))
    ans <- ans + fx * w[j + 1]
  }

  ans <- ans * h * al
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat(sprintf("       RESULTADO FINAL   I = %.6f\n", ans))
    cat("--------------------------------------------------\n")
  }
  invisible(ans)
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 4-2: FORMULAS CERRADAS DE NEWTON-COTES ===\n\n")

check_weights()

# ── Salida del libro (p. 152): f(x) = sin(x), k = 6 datos, [0, 2] ─────────
cat("--- Salida del libro: f(x) = sin(x) en [0, 2], k = 6 datos ---")
res <- newton_cotes(sin, 0, 2, k = 6)
cat(sprintf("\n Libro: 1.416117   Exacto: 1 - cos(2) = %.6f\n\n", 1 - cos(2)))

# ── Ejemplo 4.4: longitud de arco de la cardioide ─────────────────────────
# r = 2(1 + cos t),  0 <= t <= pi
# dr/dt = -2 sin t
# L = integral de sqrt(r^2 + (dr/dt)^2) dt   ->  valor exacto = 8
arc_len <- function(t) {
  r  <- 2 * (1 + cos(t))
  dr <- -2 * sin(t)
  sqrt(r^2 + dr^2)
}

cat("--- Ejemplo 4.4: longitud de arco de r = 2(1+cos t), 0 <= t <= pi ---\n")
cat("    Valor exacto: L = 8.00000\n\n")
cat(sprintf(" %8s  %14s  %14s\n", "Orden N", "Integral L", "Error"))
for (n in 2:10) {
  l <- newton_cotes(arc_len, 0, pi, k = n + 1, verbose = FALSE)
  cat(sprintf(" %8d  %14.5f  %14.2e\n", n, l, abs(8 - l)))
}
cat(sprintf(" %8s  %14.5f\n", "exacto", 8))

# ── Comparación con la tabla del libro ────────────────────────────────────
cat("\n--- Comparacion con la tabla del ejemplo 4.4 ---\n")
libro <- c("2"=8.01823, "3"=8.00803, "4"=7.99993, "5"=7.99996, "6"=8.00000,
           "7"=8.00000, "8"=8.00000, "9"=8.00197, "10"=7.99201)
cat(sprintf(" %8s  %14s  %14s\n", "Orden N", "R (esta impl.)", "Libro"))
for (n in 2:10) {
  d <- newton_cotes(arc_len, 0, pi, k = n + 1, verbose = FALSE)
  cat(sprintf(" %8d  %14.5f  %14.5f\n", n, d, libro[as.character(n)]))
}

# ── ERRATAS DETECTADAS EN LA TABLA 4.2 DEL LIBRO ──────────────────────────
# El libro atribuye a "errores de redondeo" el crecimiento del error para
# N >= 9. Eso NO es correcto: son ERRATAS en los pesos impresos. Prueba:
#
#   N = 9  el libro imprime  5788, el valor correcto es  5778
#   N = 8  el libro imprime 10946, el valor correcto es 10496
#
# Si se usa el peso erroneo 5788 en N=9 se obtiene exactamente 8.00198,
# que es el 8.00197 que reporta la tabla del libro. Lo demostramos:
demo_errata <- function() {
  nc_raw <- function(w, q, r, n) {
    h <- (pi - 0) / n
    sum(sapply(0:n, function(j) arc_len(j * h) * w[j + 1])) * h * (q / r)
  }
  w_ok  <- c(2857, 15741, 1080, 19344, 5778, 5778, 19344, 1080, 15741, 2857)
  w_bad <- c(2857, 15741, 1080, 19344, 5788, 5788, 19344, 1080, 15741, 2857)

  cat("\n--- Demostracion: las 'anomalias' son ERRATAS, no redondeo ---\n")
  cat(sprintf(" N=9 con peso 5778 (correcto): L = %.5f   sum(w) = %d\n",
              nc_raw(w_ok, 9, 89600, 9), sum(w_ok)))
  cat(sprintf(" N=9 con peso 5788 (libro)   : L = %.5f   sum(w) = %d\n",
              nc_raw(w_bad, 9, 89600, 9), sum(w_bad)))
  cat("                        el libro reporta: L = 8.00197  <- coincide\n\n")
  cat(" La prueba definitiva es la condicion alfa*sum(w) = N:\n")
  cat(sprintf("   peso 5778 -> (9/89600)*%d = %.4f  = N  OK\n",
              sum(w_ok),  9 / 89600 * sum(w_ok)))
  cat(sprintf("   peso 5788 -> (9/89600)*%d = %.4f != N  ERRATA\n",
              sum(w_bad), 9 / 89600 * sum(w_bad)))
}
demo_errata()

cat("\n Conclusion: con los pesos correctos, Newton-Cotes cerrado NO se\n")
cat(" degrada hasta N=10 (error ~1e-13 en doble precision). El riesgo real\n")
cat(" de cancelacion existe para ordenes aun mayores, donde los pesos\n")
cat(" alternan signo y crecen sin control.\n")
