# Program 6-2: Inversión de una Matriz
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 6
# Equivalente a CSL/F6-2.FOR
#
# Se arma el arreglo aumentado [A | I] y se aplica eliminacion de Gauss:
#   1) eliminacion hacia adelante con pivoteo parcial sobre las 2N columnas
#   2) sustitucion hacia atras aplicada a CADA una de las N columnas derechas
# Al terminar, las ultimas N columnas contienen A^-1.
# De paso se obtiene el determinante = producto de la diagonal (con signo).

# ── Épsilon de la máquina (igual que el listado) ──────────────────────────
# Se divide entre 2 hasta que 1 + eps ya no se distingue de 1.
machine_eps <- function() {
  eps <- 1.0
  while (1.0 + eps > 1.0) eps <- eps / 2
  eps * 2
}

# ── Subrutina GAUSS: inversión por eliminación de Gauss ───────────────────
# Recibe la matriz cuadrada A y devuelve la inversa, el determinante,
# el numero de intercambios de filas y el arreglo aumentado final.
gauss_inverse <- function(A, verbose = TRUE) {
  N  <- nrow(A)
  M  <- cbind(A, diag(N))          # arreglo aumentado [A | I], N x 2N

  eps  <- machine_eps()
  eps2 <- eps * 2
  if (verbose) cat(sprintf(" EPSILON DE LA MAQUINA = %.7e\n", eps))

  det    <- 1.0
  swaps  <- 0

  # --- Eliminacion hacia adelante con pivoteo parcial ---
  for (i in seq_len(N - 1)) {
    pv <- i
    for (j in (i + 1):N)
      if (abs(M[pv, i]) < abs(M[j, i])) pv <- j      # mayor |pivote|
    if (pv != i) {
      M[c(i, pv), ] <- M[c(pv, i), ]                 # intercambio de filas
      det   <- -det                                  # cada cambio invierte el signo
      swaps <- swaps + 1
    }
    if (M[i, i] == 0) stop("LA MATRIZ ES SINGULAR")

    # eliminacion por debajo de la diagonal (en las 2N columnas)
    for (jr in (i + 1):N) {
      if (M[jr, i] == 0) next
      r <- M[jr, i] / M[i, i]
      for (kc in (i + 1):(2 * N)) {
        temp <- M[jr, kc]
        M[jr, kc] <- M[jr, kc] - r * M[i, kc]
        # residuo de redondeo -> se fuerza a cero (como en el FORTRAN)
        if (abs(M[jr, kc]) < eps2 * abs(temp)) M[jr, kc] <- 0
      }
      M[jr, i] <- 0                                  # ya eliminado
    }
  }

  det <- det * prod(diag(M[, 1:N]))                  # producto de los pivotes
  if (M[N, N] == 0) stop("LA MATRIZ ES SINGULAR")

  # --- Sustitucion hacia atras, una vez por cada columna de la derecha ---
  for (m in (N + 1):(2 * N)) {
    M[N, m] <- M[N, m] / M[N, N]
    if (N >= 2) {
      for (nv in (N - 1):1) {
        va <- M[nv, m] - sum(M[nv, (nv + 1):N] * M[(nv + 1):N, m])
        M[nv, m] <- va / M[nv, nv]
      }
    }
  }

  list(inv = M[, (N + 1):(2 * N), drop = FALSE], det = det,
       swaps = swaps, eps = eps, aug = M)
}

# ── Impresión en el formato del libro (1P5E12.5) ──────────────────────────
print_mat <- function(X) {
  for (i in seq_len(nrow(X)))
    cat(paste(sprintf("%12.5E", X[i, ] + 0), collapse = ""), "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 6-2: INVERSION DE UNA MATRIZ ===\n")

# ── Ejemplo del libro ─────────────────────────────────────────────────────
A <- matrix(c( 2, 1, -3,
              -1, 3,  2,
               3, 1, -3), nrow = 3, byrow = TRUE)

cat("\n MATRIZ ORIGINAL\n")
print_mat(A)
cat("\n")
res <- gauss_inverse(A)
cat(sprintf("\n DETERMINANTE = %.5f\n", res$det))
cat(sprintf(" (intercambios de filas realizados: %d)\n", res$swaps))
cat("\n MATRIZ INVERSA\n")
print_mat(res$inv)

cat("\n Esperado (libro):\n")
cat("  -1.00000E+00 ~0           1.00000E+00\n")
cat("   2.72727E-01 2.72727E-01 -9.09091E-02\n")
cat("  -9.09091E-01 9.09091E-02  6.36364E-01\n")

# ── Inversa en forma racional: A^-1 = adj(A) / det(A) ─────────────────────
cat("\n Multiplicando por det(A) = 11 se ven los enteros de la adjunta:\n")
print(round(res$inv * 11))

# ── Verificación: A * A^-1 = I ────────────────────────────────────────────
cat("\n### VERIFICACION ###\n")
cat("\n A %*% A^-1 =\n")
print(round(A %*% res$inv, 12))
cat("\n A^-1 %*% A =\n")
print(round(res$inv %*% A, 12))
cat(sprintf("\n max |A A^-1 - I| = %.2e\n", max(abs(A %*% res$inv - diag(3)))))
cat(sprintf(" Comparacion con solve() de R: max diferencia = %.2e\n",
            max(abs(res$inv - solve(A)))))

# ── Uso de la inversa: resolver A x = y ───────────────────────────────────
cat("\n### USO: RESOLVER A x = y CON LA INVERSA ###\n")
y <- c(-1, 12, 0)
x <- res$inv %*% y
cat(" y = (-1, 12, 0)\n")
cat(sprintf(" x = A^-1 y = (%.6f, %.6f, %.6f)\n", x[1], x[2], x[3]))
cat(" Esperado (ejemplo de Gauss del tema): x1 = 1, x2 = 3, x3 = 2\n")

# ── Matriz singular: el programa debe detectarla ──────────────────────────
cat("\n### MATRIZ SINGULAR ###\n")
S <- matrix(c(1, 2, 3,
              2, 4, 6,
              1, 0, 1), nrow = 3, byrow = TRUE)
cat(" S (fila 2 = 2 * fila 1):\n")
print_mat(S)
tryCatch(gauss_inverse(S, verbose = FALSE),
         error = function(e) cat(" ->", conditionMessage(e), "\n"))

# ══════════════════════════════════════════════════════════════════════════════
# PROBLEMA MAL CONDICIONADO: MATRIZ DE HILBERT 4x4  (ejercicio propuesto)
# ══════════════════════════════════════════════════════════════════════════════
#   a) A^-1 (A^-1)^-1   -> deberia ser I
#   b) det(A) * det(A^-1) -> deberia ser 1
cat("\n### MATRIZ DE HILBERT (mal condicionada) ###\n")
hilbert <- function(n) outer(1:n, 1:n, function(i, j) 1 / (i + j - 1))

H     <- hilbert(4)
rH    <- gauss_inverse(H, verbose = FALSE)
Hi    <- rH$inv
rHi   <- gauss_inverse(Hi, verbose = FALSE)
Hii   <- rHi$inv                      # (A^-1)^-1, deberia volver a ser A

cat("\n H^-1 (todos enteros en teoria):\n")
print(round(Hi, 6))
cat(sprintf("\n det(H)            = %.6e\n", rH$det))
cat(sprintf(" det(H^-1)         = %.6e\n", rHi$det))
cat(sprintf(" b) det(H)*det(H^-1) = %.15f   (ideal: 1)\n", rH$det * rHi$det))

P <- Hi %*% Hii
cat("\n a) H^-1 (H^-1)^-1 =\n")
print(P)
cat(sprintf("\n max |H^-1 (H^-1)^-1 - I| = %.2e\n", max(abs(P - diag(4)))))
cat(sprintf(" max |H H^-1 - I|         = %.2e\n", max(abs(H %*% Hi - diag(4)))))
cat(sprintf(" max |(H^-1)^-1 - H|      = %.2e\n", max(abs(Hii - H))))

# Que tan rapido empeora al crecer n
cat("\n Crecimiento del error con el orden n:\n")
cat(sprintf(" %3s  %14s  %16s  %16s\n", "n", "kappa(H)", "max|H H^-1 - I|", "det*det(inv) - 1"))
for (n in c(3, 4, 6, 8, 10, 12)) {
  Hn  <- hilbert(n)
  rn  <- gauss_inverse(Hn, verbose = FALSE)
  rni <- gauss_inverse(rn$inv, verbose = FALSE)
  cat(sprintf(" %3d  %14.3e  %16.3e  %16.3e\n", n, kappa(Hn, exact = TRUE),
              max(abs(Hn %*% rn$inv - diag(n))), rn$det * rni$det - 1))
}
cat("\n Moraleja: la matriz SI tiene inversa, pero el redondeo la destruye.\n")
cat(" Cada orden de magnitud de kappa(A) cuesta ~1 cifra significativa.\n")
