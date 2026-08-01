# Program 3-2: Búsqueda de Raíces
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 3
# Equivalente a CSL/F3-2.FOR

# ── Función del ejemplo 3.2 del libro ─────────────────────────────────────
# ── -19(x-1/2)*(x-1) + e^x - e^{-2x}
# las raices están entre (0 ; 0.5) , (1 ; 1.5) , (6 ; 6.5)
func <- function(x) {
  -19 * (x - 0.5) * (x - 1) + exp(x) - exp(-2 * x)
}

# ── Búsqueda de intervalos con cambio de signo ────────────────────────────
# Recorre [a, b] con paso h e imprime los intervalos donde f cambia de signo.
root_search <- function(f, a, b, h) {
  cat(sprintf("\n INICIAL X = %.6f\n", a))
  cat(sprintf(" FINAL   X = %.6f\n", b))
  cat(sprintf(" INCREMENTO = %.6f\n\n", h))

  x  <- a
  yb <- 0.0
  i  <- 0L

  while (TRUE) {
    i <- i + 1L
    if (x > b) break
    y <- f(x)
    if (i > 1L && y * yb < 0) {
      cat(sprintf(" UN INTERVALO QUE PUEDE CONTENER UNA RAIZ: [%10.6f, %10.6f]\n",
                  x - h, x))
    }
    yb <- y
    x  <- x + h
  }

  cat("\n                *** FIN DE LA BUSQUEDA\n")
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 3-2: BUSQUEDA DE RAICES ===\n")
cat("    f(x) = -19(x-0.5)(x-1) + e^x - e^(-2x)\n")

# Caso 1: búsqueda amplia — igual al ejemplo del libro
root_search(func, a = -10, b = 10, h = 1)

# Caso 2: refinamiento en [0, 1] — igual al ejemplo del libro
root_search(func, a = 0, b = 1, h = 0.001)
