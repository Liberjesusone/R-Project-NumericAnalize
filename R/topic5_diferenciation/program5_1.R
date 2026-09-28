# Program 5-1: Cálculo de Aproximaciones por Diferencias
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 5
# Equivalente a CSL/F5-1.FOR
#
# El programa NO calcula una derivada: GENERA la formula de diferencias
# para cualquier conjunto de puntos de la reticula y cualquier orden p.
#
#   f^(p)(0) = [ a_alfa*f_alfa + ... + a_lambda*f_lambda ] / h^p  +  E
#   E = c1 * h^(L-p) * f^(L)  +  c2 * h^(L-p+1) * f^(L+1)

# ── Eliminación gaussiana con pivoteo parcial ─────────────────────────────
# Equivale a la SUBROUTINE GAUSS del listado (Cap. 6 del libro).
# Se implementa desde cero en vez de usar solve() para ser fiel al programa.
gauss_solve <- function(A, b) {
  n  <- nrow(A)
  M  <- cbind(A, b)          # matriz aumentada [A | b]

  # --- Eliminación hacia adelante ---
  for (i in seq_len(n - 1)) {
    # pivoteo parcial: buscar el mayor |elemento| en la columna i
    ipv <- i
    for (j in (i + 1):n)
      if (abs(M[ipv, i]) < abs(M[j, i])) ipv <- j
    if (ipv != i) M[c(i, ipv), ] <- M[c(ipv, i), ]   # intercambio de filas

    if (M[i, i] == 0) next
    for (jr in (i + 1):n) {
      if (M[jr, i] == 0) next
      r <- M[jr, i] / M[i, i]
      M[jr, i:(n + 1)] <- M[jr, i:(n + 1)] - r * M[i, i:(n + 1)]
    }
  }

  if (M[n, n] == 0) stop("LA MATRIZ ES SINGULAR")

  # --- Sustitución hacia atrás ---
  x <- numeric(n)
  x[n] <- M[n, n + 1] / M[n, n]
  if (n >= 2) {
    for (nv in (n - 1):1)
      x[nv] <- (M[nv, n + 1] - sum(M[nv, (nv + 1):n] * x[(nv + 1):n])) / M[nv, nv]
  }
  x
}

# ── Utilidades para expresar los coeficientes como fracciones ─────────────
# Fraccion continua: convierte un decimal en num/den exacto.
to_frac <- function(x, max_den = 1e6, tol = 1e-10) {
  if (abs(x) < tol) return(c(0, 1))
  s <- sign(x); x <- abs(x)
  h0 <- 0; h1 <- 1; k0 <- 1; k1 <- 0; b <- x
  repeat {
    a  <- floor(b)
    h2 <- a * h1 + h0;  k2 <- a * k1 + k0
    h0 <- h1; h1 <- h2; k0 <- k1; k1 <- k2
    if (k1 > max_den || abs(x - h1 / k1) < tol * max(1, x)) break
    if (abs(b - a) < tol) break
    b <- 1 / (b - a)
  }
  c(s * h1, k1)
}

gcd2 <- function(a, b) { a <- abs(a); b <- abs(b); while (b > 0) { t <- b; b <- a %% b; a <- t }; a }
lcm2 <- function(a, b) if (a == 0 || b == 0) 0 else a / gcd2(a, b) * b

# ── Núcleo: construir y resolver el sistema ───────────────────────────────
# el : vector de indices de la reticula (alfa, beta, ..., lambda)
# p  : orden de la derivada
diff_scheme <- function(el, p) {
  L <- length(el)
  if (L < p + 1) stop("Se necesitan al menos p+1 puntos")

  # Matriz de Vandermonde: fila k (potencia k-1), columna l (punto el[l])
  #   A[k, l] = el[l]^(k-1)
  # Se construyen L+2 filas: las primeras L para el sistema, las 2 ultimas
  # para calcular el termino del error.
  B <- outer(0:(L + 1), el, function(k, e) e^k)
  B[1, ] <- 1                       # fila k=1 -> el^0 = 1 (evita 0^0)

  A_sys <- B[1:L, , drop = FALSE]

  # Termino no homogeneo: cero salvo en la fila p+1, donde vale p!
  rhs <- numeric(L)
  rhs[p + 1] <- factorial(p)

  a <- gauss_solve(A_sys, rhs)      # coeficientes de la formula

  # --- Termino del error: filas L+1 y L+2 ---
  C1 <- sum(B[L + 1, ] * a)
  C2 <- sum(B[L + 2, ] * a)
  if (abs(C1) < 1e-8) C1 <- 0
  if (abs(C2) < 1e-8) C2 <- 0

  c1 <- -C1 / factorial(L)          # coeficiente de h^(L-p) f^(L)
  c2 <- -C2 / factorial(L + 1)      # coeficiente de h^(L-p+1) f^(L+1)

  list(el = el, p = p, L = L, a = a,
       c1 = c1, c2 = c2,
       C1 = C1, C2 = C2,
       pow1 = L - p, ord1 = L,
       pow2 = L - p + 1, ord2 = L + 1)
}

# ── Impresión en el formato del libro ─────────────────────────────────────
print_scheme <- function(s) {
  # Normalizacion del libro: dividir por el menor |a| no nulo
  nz <- abs(s$a[abs(s$a) > 1e-4])
  F  <- if (length(nz)) min(nz) else 1
  CF <- s$a / F

  cat("\n ESQUEMA DE DIFERENCIAS\n")
  for (k in seq_len(s$L))
    cat(sprintf(" +[%10.5f/(%8.5f H**%d)] F(%6.3fH)\n",
                CF[k], 1 / F, s$p, s$el[k]))

  cat("\n TERMINO DEL ERROR\n")
  cat(sprintf("     (%10.5f/%10.5f)H**%d  F^(%d)\n",
              -s$C1 + 0, factorial(s$L), s$pow1, s$ord1))     # +0 evita "-0.00000"
  cat(sprintf("     (%10.5f/%10.5f)H**%d  F^(%d)\n",
              -s$C2 + 0, factorial(s$L + 1), s$pow2, s$ord2))
  if (s$C1 == 0)
    cat("     (el primer termino se anula -> el error es del orden siguiente)\n")
  cat("\n--------------------------------------------------\n")
}

# ── Impresión en forma racional exacta ────────────────────────────────────
print_fraction <- function(s) {
  fr  <- lapply(s$a, to_frac)
  den <- Reduce(lcm2, vapply(fr, `[`, numeric(1), 2))
  num <- round(s$a * den)

  cat(sprintf("\n Forma racional:  f^(%d)(0) = [ ", s$p))
  cat(paste(sprintf("%+d*f(%gh)", num, s$el), collapse = " "))
  cat(sprintf(" ] / (%g h^%d)\n", den, s$p))

  fc <- to_frac(if (s$c1 != 0) s$c1 else s$c2)
  po <- if (s$c1 != 0) s$pow1 else s$pow2
  od <- if (s$c1 != 0) s$ord1 else s$ord2
  cat(sprintf(" Error dominante: E = (%d/%d) h^%d f^(%d)   ->  O(h^%d)\n",
              fc[1], fc[2], po, od, po))
}

# ── Verificación numérica de la fórmula generada ──────────────────────────
# Aplica el esquema a una funcion conocida y compara con la derivada exacta.
verify <- function(s, f, dp_exact, x0 = 0, hs = c(0.1, 0.05, 0.025)) {
  cat(sprintf("\n Verificacion con f(x) = exp(x):  f^(%d)(%g) = %.6f\n",
              s$p, x0, dp_exact))
  cat(sprintf(" %8s  %14s  %12s  %8s\n", "h", "aproximado", "error", "ratio"))
  prev <- NA
  for (h in hs) {
    ap  <- sum(s$a * f(x0 + s$el * h)) / h^s$p
    err <- abs(dp_exact - ap)
    cat(sprintf(" %8.4f  %14.8f  %12.2e  %8s\n", h, ap, err,
                if (is.na(prev)) "-" else sprintf("%.1f", prev / err)))
    prev <- err
  }
  cat(sprintf(" Orden teorico O(h^%d) -> al partir h a la mitad el error cae ~%d veces\n",
              if (s$c1 != 0) s$pow1 else s$pow2,
              2^(if (s$c1 != 0) s$pow1 else s$pow2)))
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 5-1: CALCULO DE APROXIMACIONES POR DIFERENCIAS ===\n")

# ── Ejemplo 1 del libro: 3 puntos (0,1,2), primera derivada ───────────────
cat("\n### EJEMPLO 1 DEL LIBRO ###")
cat("\n Puntos: 0, 1, 2   Orden de la derivada: 1\n")
s1 <- diff_scheme(el = c(0, 1, 2), p = 1)
print_scheme(s1)
print_fraction(s1)
cat("\n Esperado (libro): -3/2h, 4/2h, -1/2h  con E = (2/6)h^2 f'''\n")
verify(s1, exp, 1)

# ── Ejemplo 2 del libro: 5 puntos (-2..2), segunda derivada ───────────────
cat("\n\n### EJEMPLO 2 DEL LIBRO ###")
cat("\n Puntos: -2, -1, 0, 1, 2   Orden de la derivada: 2\n")
s2 <- diff_scheme(el = c(-2, -1, 0, 1, 2), p = 2)
print_scheme(s2)
print_fraction(s2)
cat("\n Esperado (libro): -1,16,-30,16,-1 sobre 12h^2  con E = (8/720)h^4 f^(6)\n")
verify(s2, exp, 1)

# ── Fórmulas clásicas reproducidas por el mismo algoritmo ─────────────────
cat("\n\n### FORMULAS CLASICAS GENERADAS POR EL MISMO ALGORITMO ###\n")

casos <- list(
  list(el = c(0, 1),           p = 1, nom = "Diferencia hacia ADELANTE (2 pts)"),
  list(el = c(-1, 0),          p = 1, nom = "Diferencia hacia ATRAS (2 pts)"),
  list(el = c(-1, 1),          p = 1, nom = "Diferencia CENTRAL (2 pts)"),
  list(el = c(0, 1, 2),        p = 1, nom = "Adelante 3 puntos"),
  list(el = c(-2, -1, 0),      p = 1, nom = "Atras 3 puntos"),
  list(el = c(-1, 0, 1),       p = 2, nom = "Central para f'' (3 pts)"),
  list(el = c(-2, -1, 0),      p = 2, nom = "Atras para f'' (3 pts)"),
  list(el = c(-2, -1, 0, 1, 2), p = 1, nom = "Central 5 puntos para f'")
)

for (cs in casos) {
  s   <- diff_scheme(cs$el, cs$p)
  fr  <- lapply(s$a, to_frac)
  den <- Reduce(lcm2, vapply(fr, `[`, numeric(1), 2))
  num <- round(s$a * den)
  po  <- if (s$c1 != 0) s$pow1 else s$pow2
  cat(sprintf("\n %s\n", cs$nom))
  cat(sprintf("   f^(%d) = [ %s ] / (%g h^%d)   O(h^%d)\n",
              s$p, paste(sprintf("%+d f(%gh)", num, s$el), collapse = " "),
              den, s$p, po))
}

# ── Retícula NO uniforme (el algoritmo también funciona) ──────────────────
cat("\n\n### RETICULA NO UNIFORME (indices decimales) ###")
cat("\n Ejemplo del enunciado: f''(0) con puntos en -2, 0.5, 1.5\n")
s3 <- diff_scheme(el = c(-2, 0.5, 1.5), p = 2)
print_scheme(s3)
print_fraction(s3)
cat("\n Nota: el algoritmo no exige espaciamiento uniforme ni indices enteros.\n")

# ══════════════════════════════════════════════════════════════════════════════
# EL PELIGRO DE h DEMASIADO PEQUEÑO
# ══════════════════════════════════════════════════════════════════════════════
# En INTEGRACION achicar h siempre mejora. En DIFERENCIACION no:
#   - error de truncamiento  ~ h      (baja al achicar h)
#   - error de redondeo      ~ eps/h  (SUBE al achicar h)
# Existe un h optimo; por debajo de el, el error CRECE.

cat("\n\n### EL PELIGRO DE h DEMASIADO PEQUENO ###\n")
cat(" f(x) = exp(x),  f'(1) = e = 2.718281828459045\n\n")
x0 <- 1
exact <- exp(1)
hs <- 10^(-(1:16))
ef <- abs((exp(x0 + hs) - exp(x0)) / hs - exact)          # adelante
ec <- abs((exp(x0 + hs) - exp(x0 - hs)) / (2 * hs) - exact)  # central
k_f <- which.min(ef)
k_c <- which.min(ec)

cat(sprintf(" %10s  %14s  %14s\n", "h", "err ADELANTE", "err CENTRAL"))
for (k in seq_along(hs)) {
  marca <- paste0(if (k == k_f) "  <- min adelante" else "",
                  if (k == k_c) "  <- min central"  else "")
  cat(sprintf(" %10.0e  %14.3e  %14.3e%s\n", hs[k], ef[k], ec[k], marca))
}

h_f <- hs[k_f]
h_c <- hs[k_c]
eps <- .Machine$double.eps
cat(sprintf("\n Epsilon de maquina: %.3e\n", eps))
cat(sprintf(" ADELANTE  O(h)  : mejor h = %.0e   teorico sqrt(eps)   = %.1e\n",
            h_f, sqrt(eps)))
cat(sprintf(" CENTRAL   O(h^2): mejor h = %.0e   teorico eps^(1/3)   = %.1e\n",
            h_c, eps^(1/3)))
cat("\n Con h = 1e-16 el error es ~2.7: el 100% del valor buscado.\n")
cat(" Moraleja: en diferenciacion, h muy pequeno es TAN malo como h grande.\n")
cat(" La integracion promedia (suaviza el ruido);\n")
cat(" la diferenciacion resta (lo amplifica).\n")
