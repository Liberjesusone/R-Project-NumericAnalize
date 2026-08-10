# Program 4-1: Reglas extendidas del Trapecio y de Simpson
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 4
# Equivalente a CSL/F4-1.FOR

# ── Función a integrar (ejemplo 4.1 del libro) ────────────────────────────
# Volumen del cuerpo de revolución de y = 1 + (x/2)^2 girado sobre el eje x
# f(x) = pi * (1 + (x/2)^2)^2       Valor exacto en [0,2]: I = 11.7286
func <- function(x) {
  pi * (1 + (x / 2)^2)^2
}

# ── Regla extendida del TRAPECIO ──────────────────────────────────────────
#   I = (h/2) * [ f0 + 2*f1 + 2*f2 + ... + 2*f(N-1) + fN ]
# Pesos: 1 en los extremos, 2 en todos los puntos interiores.
trapz <- function(f, a, b, n) {
  h <- (b - a) / n
  s <- 0
  for (i in 0:n) {
    x <- a + i * h
    w <- if (i == 0 || i == n) 1 else 2
    s <- s + w * f(x)
  }
  s * h / 2
}

# ── Regla extendida de SIMPSON ────────────────────────────────────────────
#   N par  -> solo regla 1/3:  (h/3)*[f0 + 4f1 + 2f2 + 4f3 + ... + fN]
#   N impar-> regla 3/8 en los 3 primeros intervalos + 1/3 en el resto
simps <- function(f, a, b, n) {
  h  <- (b - a) / n
  ss <- 0
  ls <- 0                      # cuántos intervalos consumió la regla 3/8

  if (n %% 2 != 0) {           # N impar: arrancar con 3/8 de Simpson
    ls <- 3
    for (i in 0:3) {
      x <- a + h * i
      w <- if (i == 0 || i == 3) 1 else 3
      ss <- ss + w * f(x)
    }
    ss <- ss * h * 3 / 8
    if (n == 3) return(ss)     # con N=3 la regla 3/8 cubre todo el dominio
  }

  # Regla 1/3 de Simpson sobre los intervalos restantes
  s <- 0
  m <- n - ls                  # número de intervalos que faltan
  for (i in 0:m) {
    x <- a + h * (i + ls)
    w <- if (i %% 2 == 1) 4 else 2
    if (i == 0 || i == m) w <- 1
    s <- s + w * f(x)
  }
  ss + s * h / 3
}

# ── Integración de Romberg (extrapolación de Richardson) ──────────────────
#   I ≈ I_h + (1/3)(I_h - I_2h)     eleva el orden de h^2 a h^4
romberg <- function(f, a, b, n) {
  i_h  <- trapz(f, a, b, n)
  i_2h <- trapz(f, a, b, n / 2)
  i_h + (i_h - i_2h) / 3
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 4-1: REGLAS DEL TRAPECIO Y DE SIMPSON ===\n")
cat("    f(x) = pi*(1 + (x/2)^2)^2   sobre [0, 2]\n")
cat("    Valor exacto: I = 11.7286\n\n")

a <- 0; b <- 2
exacto <- 11.7286

# ── Ejemplo 4.1: regla del trapecio ───────────────────────────────────────
cat("--- Ejemplo 4.1: Regla extendida del TRAPECIO ---\n")
cat(sprintf(" %5s  %10s  %12s  %12s\n", "N", "h", "I_h", "E_h"))
for (n in c(2, 4, 8, 16, 32, 64, 128)) {
  ih <- trapz(func, a, b, n)
  cat(sprintf(" %5d  %10.6f  %12.4f  %12.4f\n", n, (b - a) / n, ih, exacto - ih))
}

# ── Ejemplo 4.3: regla 1/3 de Simpson ─────────────────────────────────────
cat("\n--- Ejemplo 4.3: Regla extendida de SIMPSON (1/3) ---\n")
cat(sprintf(" %5s  %10s  %12s  %12s\n", "N", "h", "I_h", "E_h"))
for (n in c(2, 4, 8, 16, 32, 64)) {
  ih <- simps(func, a, b, n)
  cat(sprintf(" %5d  %10.6f  %12.4f  %12.4f\n", n, (b - a) / n, ih, exacto - ih))
}

# ── Salida del libro (p. 143) ─────────────────────────────────────────────
cat("\n--- Salida del libro ---\n")
cat(sprintf(" Trapecio, N = 10 : %12.5f   (libro: 11.77047)\n", trapz(func, a, b, 10)))
cat(sprintf(" Simpson,  N =  5 : %12.5f   (libro: 11.73095)\n", simps(func, a, b, 5)))
cat("   N = 5 es impar -> usa 3/8 en [0, 1.2] y 1/3 en [1.2, 2]\n")

# ── Ejemplo 4.2: integración de Romberg ───────────────────────────────────
cat("\n--- Ejemplo 4.2: Integracion de ROMBERG ---\n")
i_05  <- trapz(func, a, b, 4)    # h = 0.5
i_025 <- trapz(func, a, b, 8)    # h = 0.25
cat(sprintf(" I(h=0.50) = %.4f\n", i_05))
cat(sprintf(" I(h=0.25) = %.4f\n", i_025))
cat(sprintf(" E(h=0.25) = (1/3)(%.4f - %.4f) = %.4f\n",
            i_025, i_05, (i_025 - i_05) / 3))
cat(sprintf(" I corregida = %.4f    (equivale a trapecio con N=128)\n",
            romberg(func, a, b, 8)))

# ── Comparación de convergencia ───────────────────────────────────────────
cat("\n--- Comparacion: mismo N, distinta precision ---\n")
cat(sprintf(" %5s  %14s  %14s  %14s\n", "N", "Trapecio", "Simpson", "Romberg"))
for (n in c(4, 8, 16, 32)) {
  cat(sprintf(" %5d  %14.6f  %14.6f  %14.6f\n",
              n, trapz(func, a, b, n), simps(func, a, b, n), romberg(func, a, b, n)))
}
cat(sprintf(" %5s  %14.6f  %14.6f  %14.6f\n", "exac", exacto, exacto, exacto))
