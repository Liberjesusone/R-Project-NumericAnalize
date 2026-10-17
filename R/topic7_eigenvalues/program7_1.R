# Program 7-1: Método de Interpolación
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 7
# Equivalente a CSL/F7-1.FOR
#
# Obtiene los coeficientes del polinomio caracteristico
#   f(lambda) = det(A - lambda I) = c0 + c1 lambda + ... + cN lambda^N
# SIN desarrollar el determinante simbolicamente:
#   1) se evalua f numericamente en N+1 puntos lambda = 0, DE, 2DE, ..., N*DE
#   2) se arma la tabla de diferencias hacia adelante
#   3) la formula de Newton hacia adelante da el polinomio en s = lambda/DE;
#      los "coeficientes de Markov" pasan cada binomial C(s,k) a potencias de s
#   4) se cambia de variable s -> lambda dividiendo c_i entre DE^i
# Como f es un polinomio de grado N, N+1 puntos lo determinan EXACTAMENTE.

# ── Subrutina DETMNT: determinante por eliminación de Gauss con pivoteo ───
detmnt <- function(A, verbose = TRUE) {
  n  <- nrow(A)
  pc <- 0                                   # contador de pivoteos
  for (i in seq_len(n - 1)) {
    pv <- i
    for (j in (i + 1):n)
      if (abs(A[pv, i]) < abs(A[j, i])) pv <- j
    if (pv != i) {
      A[c(i, pv), ] <- A[c(pv, i), ]
      pc <- pc + 1
    }
    if (A[i, i] == 0) next                  # columna nula -> det = 0
    for (jr in (i + 1):n) {
      if (A[jr, i] != 0) {
        r <- A[jr, i] / A[i, i]
        A[jr, i:n] <- A[jr, i:n] - r * A[i, i:n]
      }
    }
  }
  det <- prod(diag(A))
  if (pc %% 2 == 1) det <- -det             # cada pivoteo cambia el signo
  if (verbose) {
    cat(sprintf(" DETERMINANTE = %.7e\n", det))
    cat(sprintf(" NUMERO DE PIVOTEOS = %d\n", pc))
  }
  det
}

# ── Núcleo: polinomio característico por interpolación ────────────────────
charpoly_interp <- function(B, DE, verbose = TRUE) {
  N  <- nrow(B)
  RA <- (0:N) * DE                          # lambda_0 ... lambda_N

  # 1) N+1 determinantes
  FF <- numeric(N + 1)
  for (jj in 0:N) {
    if (verbose) cat(sprintf("\n LAMBDA = %.7f\n", RA[jj + 1]))
    FF[jj + 1] <- detmnt(B - RA[jj + 1] * diag(N), verbose)
  }

  # 2) tabla de diferencias hacia adelante: DF[i, j] = Delta^j f_i
  #    (indices desplazados +1 porque R empieza en 1)
  DF <- matrix(NA_real_, N + 1, N + 1)
  DF[, 1] <- FF
  for (j in 1:N)
    for (i in 0:(N - j))
      DF[i + 1, j + 1] <- DF[i + 2, j] - DF[i + 1, j]

  if (verbose) {
    cat("\n TABLA DE DIFERENCIAS PARA LOS DETERMINANTES\n")
    for (i in 0:N)
      cat(sprintf("%8.4f", RA[i + 1]),
          sprintf("%11.3E", DF[i + 1, 1:(N - i + 1)]), "\n")
  }

  # 3) coeficientes de Markov: MV = coeficientes de C(s,K) = s(s-1)...(s-K+1)/K!
  #    en potencias de s. Se construyen multiplicando por (s-(K-1))/K.
  MV <- numeric(N + 1)                      # MV[L+1] = coef. de s^L
  CC <- numeric(N + 1)                      # CC[L+1] = coef. de s^L en f
  MV[2] <- 1                                # C(s,1) = s
  CC[1] <- DF[1, 1]                         # Delta^0 f_0
  CC[2] <- DF[1, 2]                         # Delta^1 f_0 * s
  if (verbose) cat("\n COEFICIENTES DE MARKOV\n")
  if (N >= 2) for (K in 2:N) {
    for (L in K:1) {
      MV[L + 1] <- (MV[L] - (K - 1) * MV[L + 1]) / K
      CC[L + 1] <- CC[L + 1] + MV[L + 1] * DF[1, K + 1]
    }
    if (verbose) cat(sprintf("%12.7f", MV[2:(K + 1)]), "\n")
  }

  # 4) de s a lambda:  s^i = lambda^i / DE^i
  CC <- CC / DE^(0:N)

  if (verbose) {
    cat("\n --- RESULTADO FINAL ---\n POTENCIA,    COEFICIENTES\n")
    for (i in 0:N) cat(sprintf("  %3d     %10.4f\n", i, CC[i + 1]))
  }
  list(coef = CC, FF = FF, DF = DF, RA = RA)
}

# ── Coeficientes exactos a partir de los valores propios (para comparar) ──
# det(A - lambda I) = prod(lambda_i - lambda)
poly_from_roots <- function(r) {
  p <- 1                                    # coeficientes en orden creciente
  for (ri in r) p <- c(ri * p, 0) - c(0, p) # multiplicar por (ri - lambda)
  Re(p)
}

# ── Del polinomio a las raíces: matriz de Frobenius + iteración QR ────────
# (sección 4 del tema) Para a0 + a1 x + ... + aN x^N:
#   primera fila = -a_{N-1}/a_N, ..., -a_0/a_N ; subdiagonal de unos
frobenius <- function(a) {
  N <- length(a) - 1
  p <- a / a[N + 1]
  M <- matrix(0, N, N)
  M[1, ] <- -rev(p[1:N])
  if (N >= 2) M[cbind(2:N, 1:(N - 1))] <- 1
  M
}

# Iteracion QR sin desplazamiento: A_{k+1} = R_k Q_k.
# Converge a triangular superior POR BLOQUES: bloques 1x1 (valores reales)
# y 2x2 (pares complejos conjugados).
qr_eigen <- function(M, iters = 2000, tol = 1e-9) {
  for (k in seq_len(iters)) {
    q <- qr(M)
    M <- qr.R(q) %*% qr.Q(q)
  }
  n <- nrow(M); ev <- complex(0); i <- 1
  while (i <= n) {
    if (i == n || abs(M[i + 1, i]) < tol * (abs(M[i, i]) + abs(M[i + 1, i + 1]))) {
      ev <- c(ev, complex(real = M[i, i]))
      i <- i + 1
    } else {                                # bloque 2x2
      tr <- M[i, i] + M[i + 1, i + 1]
      dt <- M[i, i] * M[i + 1, i + 1] - M[i, i + 1] * M[i + 1, i]
      ds <- sqrt(as.complex(tr^2 / 4 - dt))
      ev <- c(ev, tr / 2 + ds, tr / 2 - ds)
      i <- i + 2
    }
  }
  list(values = ev, T = M)
}

fmt_c <- function(z) {
  ifelse(abs(Im(z)) < 1e-10, sprintf("%10.6f", Re(z)),
         sprintf("%8.4f %s %.4fi", Re(z), ifelse(Im(z) < 0, "-", "+"), abs(Im(z))))
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 7-1: METODO DE INTERPOLACION ===\n")

B <- matrix(c(4, 3, 2, 1,
              3, 3, 2, 1,
              2, 2, 2, 1,
              1, 1, 1, 1), nrow = 4, byrow = TRUE)
cat("\n ORDEN DE LA MATRIZ = 4\n MATRIZ\n")
for (i in 1:4) cat(sprintf("%12.5E", B[i, ]), "\n")
cat("\n DELTA LAMBDA ? 0.5\n")

res <- charpoly_interp(B, DE = 0.5)

cat("\n Esperado (libro): 1, -7, 15, -10, 1\n")
cat("   f(lambda) = lambda^4 - 10 lambda^3 + 15 lambda^2 - 7 lambda + 1\n")

# ── Verificación independiente ────────────────────────────────────────────
cat("\n### VERIFICACION ###\n")
ev    <- eigen(B, symmetric = TRUE)$values
exact <- poly_from_roots(ev)
cat(" Coeficientes a partir de eigen():", sprintf("%.10f", exact), "\n")
cat(sprintf(" max |interpolacion - exacto| = %.2e\n", max(abs(res$coef - exact))))
cat(sprintf(" Chequeos: c0 = det(A) = %.6f   |  -c3 = traza(A) = %.6f\n",
            det(B), -res$coef[4]))

# ── Del polinomio a los valores propios (Frobenius + QR) ──────────────────
cat("\n### RAICES DEL POLINOMIO = VALORES PROPIOS ###\n")
M  <- frobenius(res$coef)
cat(" Matriz de Frobenius:\n"); print(M)
qe <- qr_eigen(M)
cat("\n Raices por iteracion QR:", fmt_c(sort(Re(qe$values))), "\n")
cat(" eigen(A) de R          :", sprintf("%10.6f", sort(ev)), "\n")
cat(" (lambda = 1 es raiz exacta: 1 - 10 + 15 - 7 + 1 = 0)\n")

# ── Ejemplo del tema: raíces de un polinomio de grado 5 ───────────────────
cat("\n### EJEMPLO DEL TEMA: x^5 - 0.2x^4 + 7x^3 + x^2 - 3.5x + 2 = 0 ###\n")
a5 <- c(2, -3.5, 1, 7, -0.2, 1)             # a0, a1, ..., a5
M5 <- frobenius(a5)
cat(" Matriz de Frobenius:\n"); print(M5)
q5 <- qr_eigen(M5)
cat("\n Matriz tras la iteracion QR (triangular por bloques):\n")
print(round(q5$T, 4))
cat("\n Raices:\n")
for (z in q5$values[order(Re(q5$values))]) cat("  ", fmt_c(z), "\n")
cat(" Esperado (tema): -0.9085, 0.1403 +/- 2.7314i, 0.4139 +/- 0.3506i\n")
cat(sprintf(" Comprobacion con polyroot(): max diferencia = %.2e\n",
            max(abs(sort(Mod(q5$values)) - sort(Mod(polyroot(a5)))))))

# ── Sistema homogéneo del inicio del tema ─────────────────────────────────
cat("\n### SISTEMA HOMOGENEO DEL TEMA ###\n")
H <- matrix(c(3, -2, 1,
              1,  1, 2,
              4, -1, 3), nrow = 3, byrow = TRUE)
cat(" det(A) =", sprintf("%.2e", detmnt(H, verbose = FALSE)),
    " -> hay soluciones no triviales\n")
cat(" A %*% (1, 1, -1) =", H %*% c(1, 1, -1), "  -> x = beta, y = beta, z = -beta\n")

# ══════════════════════════════════════════════════════════════════════════════
# SENSIBILIDAD AL INTERVALO DELTA LAMBDA
# ══════════════════════════════════════════════════════════════════════════════
# Delta^N f / DE^N aproxima f^(N)/... como en el capitulo 5: si DE es muy
# pequeno, las diferencias restan numeros casi iguales y se pierden cifras.
cat("\n### SENSIBILIDAD AL INTERVALO DELTA LAMBDA ###\n")
cat(sprintf(" %10s  %18s  %14s\n", "DE", "max|error coef|", "coef. lambda^4"))
for (de in c(10, 2, 0.5, 1e-1, 1e-2, 1e-3, 1e-4)) {
  r <- charpoly_interp(B, DE = de, verbose = FALSE)
  cat(sprintf(" %10.0e  %18.3e  %14.6f\n", de, max(abs(r$coef - exact)), r$coef[5]))
}
cat("\n Moraleja: igual que en diferenciacion numerica, un DE demasiado pequeno\n")
cat(" destruye la tabla de diferencias por cancelacion.\n")
