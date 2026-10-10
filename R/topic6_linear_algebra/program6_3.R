# Program 6-3: Descomposición LU con Pivoteo
# Nakamura, Metodos Numericos Aplicados con Software, Cap. 6
# Equivalente a CSL/F6-3.FOR
#
# Factoriza A = L U  (metodo de Doolittle: L con unos en la diagonal).
# L y U se guardan juntas en una sola matriz "compacta" EL:
#   - sobre y en la diagonal  -> U
#   - debajo de la diagonal   -> L  (los unos de L no se guardan)
# Con pivoteo, lo que se factoriza es P A, y el vector IP guarda la permutacion.
# det(A) = (+/-) producto de la diagonal de U.

# ── Subrutina PIVOT ───────────────────────────────────────────────────────
# Busca en la columna j (filas j..N) el mayor |A(k,j)| e intercambia filas
# tanto en A como en EL (las partes de L ya calculadas viajan con su fila).
pivot <- function(st, j) {
  N  <- nrow(st$A)
  jj <- j - 1 + which.max(abs(st$A[j:N, j]))
  if (jj != j) {
    st$A[c(j, jj), ]  <- st$A[c(jj, j), ]
    st$EL[c(j, jj), ] <- st$EL[c(jj, j), ]
    st$IP[c(j, jj)]   <- st$IP[c(jj, j)]
    st$swaps <- st$swaps + 1
    cat(sprintf(" Pivoteo en columna %d: se intercambian filas %d y %d\n", j, j, jj))
  }
  st
}

# ── Descomposición LU (Doolittle), renglón de U / columna de L ────────────
lu_decomp <- function(A, pivoting = TRUE) {
  N  <- nrow(A)
  st <- list(A = A, EL = matrix(0, N, N), IP = 1:N, swaps = 0)

  if (pivoting) st <- pivot(st, 1)

  # a) primer renglon de U
  st$EL[1, ] <- st$A[1, ]
  if (st$EL[1, 1] == 0) stop("PIVOTE CERO (u11 = 0): l21 = a21/0, active el pivoteo")
  # b) primera columna de L
  if (N >= 2) st$EL[2:N, 1] <- st$A[2:N, 1] / st$EL[1, 1]

  if (N >= 2) for (m in 2:N) {
    if (pivoting) st <- pivot(st, m)

    # e) m-esimo renglon de U:  u_mj = a_mj - sum_k l_mk u_kj
    for (j in m:N) {
      s <- sum(st$EL[m, 1:(m - 1)] * st$EL[1:(m - 1), j])
      st$EL[m, j] <- st$A[m, j] - s
    }
    if (st$EL[m, m] == 0 && m < N) stop(sprintf("PIVOTE CERO (u%d%d = 0): active el pivoteo", m, m))

    # f) m-esima columna de L:  l_im = (a_im - sum_k l_ik u_km) / u_mm
    if (m < N) for (i in (m + 1):N) {
      s <- sum(st$EL[i, 1:(m - 1)] * st$EL[1:(m - 1), m])
      st$EL[i, m] <- (st$A[i, m] - s) / st$EL[m, m]
    }
  }

  # determinante: producto de la diagonal de U, signo segun # de intercambios
  de <- prod(diag(st$EL))
  if (st$swaps %% 2 == 1) de <- -de

  L <- st$EL; L[upper.tri(L, diag = TRUE)] <- 0; diag(L) <- 1
  U <- st$EL; U[lower.tri(U)] <- 0

  list(EL = st$EL, L = L, U = U, IP = st$IP, swaps = st$swaps,
       det = de, PA = st$A)
}

# ── Resolver A x = y usando la factorización ──────────────────────────────
#   L z = P y  (sustitucion hacia adelante)
#   U x = z    (sustitucion hacia atras)
lu_solve <- function(f, y) {
  N <- nrow(f$EL)
  b <- y[f$IP]                         # aplicar la misma permutacion a y
  z <- numeric(N)
  z[1] <- b[1]
  if (N >= 2) for (i in 2:N)
    z[i] <- b[i] - sum(f$EL[i, 1:(i - 1)] * z[1:(i - 1)])
  x <- numeric(N)
  x[N] <- z[N] / f$EL[N, N]
  if (N >= 2) for (i in (N - 1):1)
    x[i] <- (z[i] - sum(f$EL[i, (i + 1):N] * x[(i + 1):N])) / f$EL[i, i]
  list(z = z, x = x)
}

print_mat <- function(X) {
  for (i in seq_len(nrow(X)))
    cat(paste(sprintf("%12.5E", X[i, ] + 0), collapse = ""), "\n")
}

show_lu <- function(f) {
  cat(" PERMUTACION\n")
  cat(sprintf("%3d", f$IP), "\n")
  cat("\n MATRICES LU EN FORMA COMPACTA\n")
  print_mat(f$EL)
  cat(sprintf("\n DETERMINANTE = %.5f\n", f$det))
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== PROGRAMA 6-3: DESCOMPOSICION LU ===\n")

A <- matrix(c( 2, 1, -3,
              -1, 3,  2,
               3, 1, -3), nrow = 3, byrow = TRUE)
cat("\n N = 3\n MATRIZ ORIGINAL\n")
print_mat(A)

# ── Caso 1: CON pivoteo (salida del libro) ────────────────────────────────
cat("\n### CON PIVOTEO (NP = 1) ###\n")
f1 <- lu_decomp(A, pivoting = TRUE)
show_lu(f1)
cat("\n Esperado (libro): permutacion 3 2 1\n")
cat("   3.00000E+00  1.00000E+00 -3.00000E+00\n")
cat("  -3.33333E-01  3.33333E+00  1.00000E+00\n")
cat("   6.66667E-01  1.00000E-01 -1.10000E+00     DETERMINANTE = 11\n")

cat("\n Separadas:\n L =\n"); print(round(f1$L, 5))
cat(" U =\n");               print(round(f1$U, 5))
cat(sprintf("\n Verificacion: max |L U - P A| = %.2e\n",
            max(abs(f1$L %*% f1$U - A[f1$IP, ]))))

# ── Caso 2: SIN pivoteo (ejemplo de la teoría) ────────────────────────────
cat("\n### SIN PIVOTEO (NP = 0) ###\n")
f0 <- lu_decomp(A, pivoting = FALSE)
show_lu(f0)
cat("\n Esperado (teoria): L = [1 0 0; -0.5 1 0; 1.5 -0.1429 1]\n")
cat("                    U = [2 1 -3; 0 3.5 0.5; 0 0 1.5714]\n")
cat(sprintf(" det = 2 * 3.5 * 1.5714 = %.5f\n", prod(diag(f0$U))))
cat(sprintf("\n Verificacion: max |L U - A| = %.2e\n",
            max(abs(f0$L %*% f0$U - A))))

# ── Uso: resolver A x = y con la factorización ────────────────────────────
cat("\n### USO: RESOLVER A x = y CON LU ###\n")
y <- c(-1, 12, 0)
cat(" y = (-1, 12, 0)\n")
for (cs in list(list(f = f0, nom = "sin pivoteo"), list(f = f1, nom = "con pivoteo"))) {
  s <- lu_solve(cs$f, y)
  cat(sprintf("\n %s:\n", cs$nom))
  cat(sprintf("   L z = P y  ->  z = (%s)\n", paste(sprintf("%.5f", s$z), collapse = ", ")))
  cat(sprintf("   U x = z    ->  x = (%s)\n", paste(sprintf("%.5f", s$x), collapse = ", ")))
}
cat("\n Esperado: x1 = 1, x2 = 3, x3 = 2 (igual que por Gauss)\n")

# Ventaja de LU: se factoriza UNA vez y se reutiliza para muchos y
cat("\n Misma factorizacion, otros lados derechos (sin refactorizar):\n")
for (yy in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
  s <- lu_solve(f1, yy)
  cat(sprintf("   y = (%s)  ->  x = (%s)\n", paste(yy, collapse = ","),
              paste(sprintf("%9.5f", s$x), collapse = " ")))
}
cat(" (son las columnas de A^-1: coincide con el Programa 6-2)\n")

# ── Por qué hace falta el pivoteo ─────────────────────────────────────────
cat("\n### POR QUE HACE FALTA EL PIVOTEO ###\n")
B <- matrix(c(0, 10,  1,
              1,  3, -1,
              2,  4,  1), nrow = 3, byrow = TRUE)
cat(" Matriz con a11 = 0 (ejemplo de pivoteo del tema):\n")
print_mat(B)
cat("\n Sin pivoteo:\n")
tryCatch(lu_decomp(B, pivoting = FALSE),
         error = function(e) cat(" ->", conditionMessage(e), "\n"))
cat("\n Con pivoteo:\n")
fb <- lu_decomp(B, pivoting = TRUE)
show_lu(fb)
sb <- lu_solve(fb, c(2, 6, 5))
cat(sprintf("\n Con y = (2, 6, 5):  x = (%s)\n",
            paste(sprintf("%.4f", sb$x), collapse = ", ")))
cat(" Esperado (teoria): x1 = 2.7187, x2 = 0.4062, x3 = -2.0625\n")
