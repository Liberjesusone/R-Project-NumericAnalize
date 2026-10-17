# Program 7-2: Householder / Bisección
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 7
# Equivalente a CSL/F7-2.FOR
#
# Valores propios de una matriz REAL SIMETRICA en dos etapas:
#   ETAPA 1: Householder -> A se reduce a tridiagonal con N-2 transformaciones
#            de similaridad A <- P A P, P = I - u u^T / h  (no cambia los
#            valores propios)
#   ETAPA 2: el determinante de una tridiagonal sale de una recurrencia de
#            3 terminos (O(N)); sus raices se buscan por biseccion + interpolacion
#            lineal, usando que los valores propios de la submatriz principal
#            L-1 SEPARAN a los de la submatriz L (entrelazamiento).

# ── ETAPA 1: reducción de Householder ─────────────────────────────────────
householder_tridiag <- function(A, verbose = TRUE) {
  N <- nrow(A)
  for (ir in seq_len(N - 2)) {
    idx <- (ir + 1):N
    s   <- sum(A[idx, ir]^2)                # S = sum a_{i,ir}^2,  i > ir
    if (s == 0) next                        # la columna ya esta reducida
    w   <- if (A[ir + 1, ir] < 0) -1 else 1 # signo(a_{ir+1,ir})
    ssr <- sqrt(s)                          # |G|
    h   <- s + abs(A[ir + 1, ir]) * ssr     # h = G^2 + G a_{ir+1,ir}

    u <- numeric(N)
    u[idx]    <- A[idx, ir]
    u[ir + 1] <- A[ir + 1, ir] + ssr * w    # a_{ir+1,ir} + G

    # P A P  sin construir P (como el FORTRAN):
    #   = A - (u (Au)^T + (Au) u^T)/h + (u^T A u / h^2) u u^T
    Au  <- as.vector(A %*% u)
    uAu <- sum(u * Au)
    A <- A - (outer(u, Au) + outer(Au, u)) / h + (uAu / h^2) * outer(u, u)
    A[abs(A) < 1e-6] <- 0                   # limpieza de residuos (umbral del libro)

    if (verbose) {
      cat(sprintf("\n Paso %d:  G = %.6f   h = %.6f   u = (%s)\n", ir, ssr * w, h,
                  paste(sprintf("%.5f", u), collapse = ", ")))
    }
  }
  A
}

# ── DETERM: determinante de la submatriz principal L x L tridiagonal ──────
# p_0 = 1,  p_1 = a_1 - x,  p_k = (a_k - x) p_{k-1} - b_{k-1}^2 p_{k-2}
tri_det <- function(a, b, x, L) {
  p0 <- 1
  p1 <- a[1] - x
  if (L == 1) return(p1)
  for (k in 2:L) {
    p2 <- (a[k] - x) * p1 - b[k - 1]^2 * p0
    p0 <- p1; p1 <- p2
  }
  p1
}

# ── BISEC: bisección mientras el intervalo es > 1, después interpolación ──
bisec <- function(a, b, L, xl, xh, tol = 1e-12, maxit = 500) {
  yl <- tri_det(a, b, xl, L)
  yh <- tri_det(a, b, xh, L)
  xm <- xl
  for (it in seq_len(maxit)) {
    xb <- xm
    if (xh - xl > 1) {
      xm <- (xl + xh) / 2                   # biseccion
    } else {
      xm <- (xl * yh - xh * yl) / (yh - yl) # interpolacion lineal (falsa posicion)
      # si el redondeo deja xm fuera de (xl, xh), se vuelve a biseccion
      if (!is.finite(xm) || xm <= xl || xm >= xh) xm <- (xl + xh) / 2
    }
    ym <- tri_det(a, b, xm, L)
    if (ym == 0 || xh - xl < tol || (xh - xl <= 1 && abs(xb - xm) < tol)) break
    if (yl * ym < 0) { xh <- xm; yh <- ym } else { xl <- xm; yl <- ym }
  }
  c(root = unname(xm), iters = it)
}

# ── ETAPA 2: valores propios de todas las submatrices principales ─────────
tridiag_eigen <- function(Tm, verbose = TRUE) {
  N <- nrow(Tm)
  a <- diag(Tm)
  b <- if (N >= 2) Tm[cbind(2:N, 1:(N - 1))] else numeric(0)
  # cota de Gershgorin: todos los valores propios estan en [-R, R]
  R <- max(abs(a) + c(0, abs(b)) + c(abs(b), 0)) + 1

  EI <- vector("list", N)
  EI[[1]] <- a[1]
  if (verbose) cat(sprintf("\n VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 1 x 1\n%12.5E\n", a[1]))
  if (N >= 2) for (L in 2:N) {
    ext <- c(-R, EI[[L - 1]], R)            # los de L-1 separan a los de L
    EI[[L]] <- sapply(1:L, function(j) unname(bisec(a, b, L, ext[j], ext[j + 1])["root"]))
    if (verbose) {
      if (L < N) cat(sprintf("\n VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE %d x %d\n", L, L))
      else cat("\n RESULTADOS FINALES (VALORES PROPIOS DE LA MATRIZ COMPLETA)\n",
               "---------------------------------------------------------\n")
      cat(sprintf("%12.5E", EI[[L]]), "\n")
    }
  }
  EI
}

print_mat <- function(X) for (i in seq_len(nrow(X))) cat(sprintf("%12.5E", X[i, ] + 0), "\n")

# ── Métodos de las potencias (sección 3 del tema) ─────────────────────────
# Directo:  v = A u,       lambda = (v,v)/(v,u),  u <- v/lambda  -> mayor |lambda|
# Inverso:  A v = u,       mu     = (v,v)/(v,u),  u <- v/mu      -> menor |lambda| = 1/mu
power_method <- function(A, u, iters, inverse = FALSE, show = integer(0)) {
  lam <- NA
  for (k in seq_len(iters)) {
    v   <- if (inverse) solve(A, u) else as.vector(A %*% u)
    mu  <- sum(v * v) / sum(v * u)
    u   <- v / mu
    lam <- if (inverse) 1 / mu else mu
    if (k %in% show)
      cat(sprintf(" %3d  %12.6f   %s\n", k, lam, paste(sprintf("%10.7f", u), collapse = " ")))
  }
  list(value = lam, vector = u)
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 7-2: HOUSEHOLDER / BISECCION ===\n")
cat(" (SOLO PARA MATRICES SIMETRICAS)\n")

A <- matrix(c(4, 3, 2, 1,
              3, 3, 2, 1,
              2, 2, 2, 1,
              1, 1, 1, 1), nrow = 4, byrow = TRUE)
cat("\n MATRIZ ORIGINAL (DEBE SER SIMETRICA)\n")
print_mat(A)

Tm <- householder_tridiag(A)
cat("\n MATRIZ DE HESSENBERG O TRIDIAGONAL\n")
print_mat(Tm)

EI <- tridiag_eigen(Tm)
cat(" ---------------------------------------------------------\n")
cat("\n Esperado (libro):\n")
cat("   1x1: 4.00000E+00\n")
cat("   2x2: 7.25082E-01  8.27492E+00\n")
cat("   3x3: 3.81085E-01  9.94722E-01  8.29086E+00\n")
cat("   4x4: 2.83250E-01  4.26021E-01  9.99996E-01  8.29086E+00\n")

# ── Verificación ──────────────────────────────────────────────────────────
cat("\n### VERIFICACION ###\n")
cat(sprintf(" Traza:  A = %.6f   T = %.6f   (se conserva)\n", sum(diag(A)), sum(diag(Tm))))
cat(sprintf(" Determinante:  A = %.6f   T = %.6f\n", det(A), det(Tm)))
cat(" eigen(A) de R:", sprintf("%12.8f", sort(eigen(A, symmetric = TRUE)$values)), "\n")
cat(" Householder+biseccion:", sprintf("%12.8f", EI[[4]]), "\n")
cat(sprintf(" max diferencia = %.2e\n", max(abs(EI[[4]] - sort(eigen(A)$values)))))

cat("\n Entrelazamiento (cada fila separa a la siguiente):\n")
for (L in 1:4) cat(sprintf("   %dx%d: ", L, L), sprintf("%9.5f", EI[[L]]), "\n")

# ── Ejemplo del tema: matriz 3x3 ──────────────────────────────────────────
cat("\n### EJEMPLO DEL TEMA: A = [3 1 4; 1 7 2; 4 2 0] ###\n")
A3 <- matrix(c(3, 1, 4,
               1, 7, 2,
               4, 2, 0), nrow = 3, byrow = TRUE)
T3 <- householder_tridiag(A3)
cat("\n Tridiagonal:\n"); print_mat(T3)
cat(" Esperado (tema): [3 -4.1231 0; -4.1231 1.3529 3.4118; 0 3.4118 5.6471]\n")

# P es simetrica y ortogonal: P P = I
u <- c(0, 1 + sqrt(17), 4); h <- 17 + sqrt(17)
P <- diag(3) - outer(u, u) / h
cat(sprintf("\n P = I - u u^T/h con u = (0, %.4f, 4), h = %.4f\n", u[2], h))
cat(sprintf(" max|P P - I| = %.2e     max|P - P^T| = %.2e\n",
            max(abs(P %*% P - diag(3))), max(abs(P - t(P)))))

E3 <- tridiag_eigen(T3, verbose = FALSE)
cat("\n p1 = 0  ->", sprintf("%9.4f", E3[[1]]), "\n")
cat(" p2 = 0  ->", sprintf("%9.4f", E3[[2]]), "   (tema: -2.0289, 6.3810)\n")
cat(" p3 = 0  ->", sprintf("%9.4f", E3[[3]]), "   (tema: -2.8941, 4.3861, 8.508)\n")

# ── Aplicación: 3 masas y 4 resortes (ejemplo físico del tema) ────────────
# Con y_j = exp(i w t) f_j  ->  A f = w^2 f,  A = K/m  (simetrica, tridiagonal)
# Todos los k = 1 y m = 1:
cat("\n### APLICACION: SISTEMA MASA-RESORTE (k = 1, m = 1) ###\n")
K <- matrix(c( 2, -1,  0,
              -1,  2, -1,
               0, -1,  2), nrow = 3, byrow = TRUE)
EK <- tridiag_eigen(householder_tridiag(K, verbose = FALSE), verbose = FALSE)[[3]]
cat(" lambda = w^2 :", sprintf("%9.5f", EK), "\n")
cat(" w (rad/s)    :", sprintf("%9.5f", sqrt(EK)), "\n")
cat(" Exacto: 2 - sqrt(2), 2, 2 + sqrt(2) =", sprintf("%.5f", c(2 - sqrt(2), 2, 2 + sqrt(2))), "\n")
cat(" Modos (vectores propios):\n")
print(round(eigen(K, symmetric = TRUE)$vectors[, 3:1], 4))
cat(" Modo 1: las 3 masas en fase | Modo 2: la del centro quieta | Modo 3: alternadas\n")

# ══════════════════════════════════════════════════════════════════════════════
# EXTRA: MÉTODO DE LAS POTENCIAS (tabla 7.3 del tema)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n### METODO DE LAS POTENCIAS (tabla 7.3) ###\n")
cat("   k        lambda    elementos de u\n")
pm <- power_method(A3, c(1, 1, 1), 20, show = c(1:5, 10, 15, 18:20))
cat(sprintf("\n Mayor valor propio: %.6f   (Householder+biseccion: %.6f)\n",
            pm$value, max(E3[[3]])))

cat("\n### POTENCIAS INVERSAS (menor |lambda|) ###\n")
cat("   k        lambda    elementos de u\n")
ip <- power_method(A3, c(1, 1, 1), 40, inverse = TRUE, show = c(1:5, 10, 20, 30, 40))
cat(sprintf("\n Menor |valor propio|: %.6f   (Householder+biseccion: %.6f)\n",
            ip$value, E3[[3]][which.min(abs(E3[[3]]))]))

cat("\n Con la matriz 4x4 del programa:\n")
p4 <- power_method(A, c(1, 1, 1, 1), 30)
i4 <- power_method(A, c(1, 1, 1, 1), 60, inverse = TRUE)
cat(sprintf("   potencias:          %.8f   (mayor)\n", p4$value))
cat(sprintf("   potencias inversas: %.8f   (menor)\n", i4$value))

cat("\n Velocidad: el error cae como |lambda_2/lambda_1|^k\n")
ev3 <- sort(abs(eigen(A3)$values), decreasing = TRUE)
ev4 <- sort(abs(eigen(A)$values),  decreasing = TRUE)
cat(sprintf("   3x3: |%.4f/%.4f| = %.3f     4x4: |%.4f/%.4f| = %.3f\n",
            ev3[2], ev3[1], ev3[2] / ev3[1], ev4[2], ev4[1], ev4[2] / ev4[1]))
