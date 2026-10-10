# Guía de Álgebra Lineal Numérica — Capítulo 6 (Nakamura)

> Universidad de Los Andes · CeSiMo · Análisis Numérico
> Programas 6-2 (Inversión) y 6-3 (Descomposición LU) implementados en R

---

# PARTE I — FORMULARIO

## 1. Eliminación de Gauss

**Eliminación hacia adelante** (para cada columna $k$, cada fila $i>k$):

$$r = \frac{a_{i,k}}{a_{k,k}}, \qquad a_{i,j} \leftarrow a_{i,j} - r\,a_{k,j}, \qquad y_i \leftarrow y_i - r\,y_k$$

**Sustitución hacia atrás:**

$$x_N = \frac{y_N}{a_{N,N}}, \qquad x_i = \frac{y_i - \sum_{j=i+1}^{N} a_{i,j}x_j}{a_{i,i}}, \quad i = N-1,\dots,1$$

Al elemento $a_{k,k}$ por el que se divide se le llama **pivote**.

---

## 2. Pivoteo parcial

Antes de eliminar la columna $k$, se busca la fila $p \ge k$ con el **mayor $|a_{p,k}|$**
y se intercambian las filas $k$ y $p$.

- Evita dividir entre cero ($a_{k,k}=0$).
- Evita dividir entre números **muy pequeños**, que amplifican el redondeo.
- Cada intercambio **cambia el signo** del determinante.

---

## 3. Inversa por Gauss-Jordan

$$\bigl[\,A \mid I\,\bigr] \;\xrightarrow{\text{operaciones de fila}}\; \bigl[\,I \mid A^{-1}\,\bigr]$$

Resolver $AX = I$ es resolver $N$ sistemas $A\mathbf{x}_j = \mathbf{e}_j$ a la vez:
la columna $j$ de $A^{-1}$ es la solución con lado derecho $\mathbf{e}_j$.

---

## 4. Descomposición LU (Doolittle)

$$A = LU, \qquad L = \begin{bmatrix}1&0&0\\l_{21}&1&0\\l_{31}&l_{32}&1\end{bmatrix},\quad
U = \begin{bmatrix}u_{11}&u_{12}&u_{13}\\0&u_{22}&u_{23}\\0&0&u_{33}\end{bmatrix}$$

Se calcula alternando **renglón de $U$ → columna de $L$**:

$$\boxed{\;u_{n,j} = a_{n,j} - \sum_{k=1}^{n-1} l_{n,k}\,u_{k,j}, \quad j = n,\dots,N\;}$$

$$\boxed{\;l_{i,n} = \frac{a_{i,n} - \sum_{k=1}^{n-1} l_{i,k}\,u_{k,n}}{u_{n,n}}, \quad i = n+1,\dots,N\;}$$

Con $n=1$ las sumas son vacías: $u_{1,j} = a_{1,j}$ y $l_{i,1} = a_{i,1}/u_{1,1}$.

**Forma compacta** (lo que imprime el programa): una sola matriz con $U$ arriba
(incluida la diagonal) y $L$ abajo; los unos de $L$ no se guardan.

$$\begin{bmatrix} u_{11} & u_{12} & u_{13} \\ l_{21} & u_{22} & u_{23} \\ l_{31} & l_{32} & u_{33} \end{bmatrix}$$

---

## 5. Resolver $A\mathbf{x} = \mathbf{y}$ con LU

$$LU\mathbf{x} = \mathbf{y} \;\Rightarrow\; \begin{cases} L\mathbf{z} = \mathbf{y} & \text{(sustitución hacia adelante)} \\ U\mathbf{x} = \mathbf{z} & \text{(sustitución hacia atrás)} \end{cases}$$

$$z_1 = y_1, \qquad z_i = y_i - \sum_{j=1}^{i-1} l_{i,j}\,z_j$$

$$x_N = \frac{z_N}{u_{N,N}}, \qquad x_i = \frac{z_i - \sum_{j=i+1}^{N} u_{i,j}\,x_j}{u_{i,i}}$$

**Con pivoteo** se factoriza $PA = LU$, así que antes hay que permutar $\mathbf{y}$:
$L\mathbf{z} = P\mathbf{y}$.

---

## 6. Determinantes

| Método | Fórmula | Costo |
|---|---|:---:|
| Definición (permutaciones) | $\sum \pm\, a_{i1}a_{j2}\cdots a_{rN}$ | $N!$ términos |
| $2\times2$ | $a_{11}a_{22} - a_{12}a_{21}$ | — |
| Sarrus ($3\times3$) | 3 diagonales $(+)$ − 3 diagonales $(-)$ | solo $3\times3$ |
| **LU / Gauss** | $\det A = (-1)^{s}\prod_{k} u_{k,k}$ | $\sim \tfrac23 N^3$ |

donde $s$ = número de intercambios de filas. Se basa en dos reglas:
$\det(BC) = \det B \det C$ y el determinante de una triangular es el producto de su diagonal.
Como $\det L = 1$, queda $\det A = \det U$.

---

## 7. Costo computacional (operaciones aproximadas)

| Tarea | Operaciones |
|---|:---:|
| Gauss (un sistema) | $\tfrac23 N^3$ |
| Factorizar LU | $\tfrac23 N^3$ |
| Resolver con LU ya factorizada | $2N^2$ por cada $\mathbf{y}$ |
| Invertir $A$ | $\approx 2N^3$ |
| Regla de Cramer / determinante por definición | $\sim N!$ |

> Con $N = 20$: $\tfrac23 N^3 \approx 5\,300$ operaciones contra $20! \approx 2.4\times10^{18}$.

---

# PARTE II — TEORÍA CON EJEMPLOS RESUELTOS

Todos los ejemplos usan la misma matriz del libro:

$$A = \begin{bmatrix} 2 & 1 & -3 \\ -1 & 3 & 2 \\ 3 & 1 & -3 \end{bmatrix}, \qquad \mathbf{y} = \begin{bmatrix} -1 \\ 12 \\ 0 \end{bmatrix}$$

## Ejemplo 1 — Eliminación de Gauss (sin pivoteo)

**Columna 1:** $r_2 = -\tfrac12$, $r_3 = \tfrac32$

$$\left[\begin{array}{ccc|c} 2 & 1 & -3 & -1 \\ 0 & 7/2 & 1/2 & 23/2 \\ 0 & -1/2 & 3/2 & 3/2 \end{array}\right]$$

**Columna 2:** $r_3 = \dfrac{-1/2}{7/2} = -\tfrac17$

$$\left[\begin{array}{ccc|c} 2 & 1 & -3 & -1 \\ 0 & 7/2 & 1/2 & 23/2 \\ 0 & 0 & 11/7 & 22/7 \end{array}\right]$$

**Sustitución hacia atrás:**

$$x_3 = \frac{22/7}{11/7} = 2, \qquad x_2 = \frac{23/2 - \tfrac12(2)}{7/2} = 3, \qquad x_1 = \frac{-1 - 3 + 3(2)}{2} = 1$$

> **Fíjate:** los multiplicadores $r$ que usamos ($-\tfrac12,\ \tfrac32,\ -\tfrac17$) son
> exactamente los $l_{ij}$ de la descomposición LU, y la matriz triangular final es $U$.
> **Gauss y LU son el mismo cálculo**: LU solo guarda los multiplicadores.

---

## Ejemplo 2 — Pivoteo

$$\left[\begin{array}{ccc|c} 0 & 10 & 1 & 2 \\ 1 & 3 & -1 & 6 \\ 2 & 4 & 1 & 5 \end{array}\right]$$

Con $a_{11} = 0$, el primer paso sería $r = 1/0$. Se intercambian filas para poner
el mayor $|a_{i1}|$ (el $2$) arriba, se elimina, y en la columna 2 se vuelve a
pivotear (el $10$ sube):

$$\left[\begin{array}{ccc|c} 2 & 4 & 1 & 5 \\ 0 & 10 & 1 & 2 \\ 0 & 0 & -1.6 & 3.3 \end{array}\right]
\;\Rightarrow\; x_3 = -2.0625,\; x_2 = 0.40625,\; x_1 = 2.71875$$

> El pivoteo no cambia la solución: intercambiar dos ecuaciones no cambia el sistema.

---

## Ejemplo 3 — Inversa por $[A \mid I]$

$$\left[\begin{array}{ccc|ccc} 2 & 1 & -3 & 1 & 0 & 0 \\ -1 & 3 & 2 & 0 & 1 & 0 \\ 3 & 1 & -3 & 0 & 0 & 1 \end{array}\right]
\;\longrightarrow\;
\left[\begin{array}{ccc|ccc} 1 & 0 & 0 & -1 & 0 & 1 \\ 0 & 1 & 0 & 0.2727 & 0.2727 & -0.0909 \\ 0 & 0 & 1 & -0.9091 & 0.0909 & 0.6364 \end{array}\right]$$

Las mismas operaciones de fila se aplican a las 6 columnas. En forma exacta:

$$A^{-1} = \frac{1}{11}\begin{bmatrix} -11 & 0 & 11 \\ 3 & 3 & -1 \\ -10 & 1 & 7 \end{bmatrix}$$

El $11$ es $\det A$: por eso los decimales son onceavos ($0.2727 = 3/11$).

**¿Por qué funciona?** Hacer operaciones de fila equivale a multiplicar por la
izquierda por una matriz $G$. Si al final $GA = I$, entonces $G = A^{-1}$, y como
la parte derecha empezó siendo $I$, terminó siendo $GI = A^{-1}$.

---

## Ejemplo 4 — LU sin pivoteo

| Paso | Cálculo | Resultado |
|---|---|---|
| a) renglón 1 de $U$ | $u_{1j} = a_{1j}$ | $2,\ 1,\ -3$ |
| b) columna 1 de $L$ | $l_{21} = -1/2,\ l_{31} = 3/2$ | $-0.5,\ 1.5$ |
| c) renglón 2 de $U$ | $u_{22} = 3 - (-0.5)(1)$, $u_{23} = 2 - (-0.5)(-3)$ | $3.5,\ 0.5$ |
| d) columna 2 de $L$ | $l_{32} = \dfrac{1 - (1.5)(1)}{3.5}$ | $-0.142857$ |
| e) renglón 3 de $U$ | $u_{33} = -3 - (1.5)(-3) - (-0.142857)(0.5)$ | $1.571429$ |

$$\text{compacta} = \begin{bmatrix} 2 & 1 & -3 \\ -0.5 & 3.5 & 0.5 \\ 1.5 & -0.1429 & 1.5714 \end{bmatrix}, \qquad \det A = 2 \cdot 3.5 \cdot 1.5714 = 11$$

**Resolviendo $A\mathbf{x} = \mathbf{y}$:**

$$L\mathbf{z} = \mathbf{y}: \quad z_1 = -1,\quad z_2 = 12 - (-0.5)(-1) = 11.5,\quad z_3 = 0 - 1.5(-1) - (-0.142857)(11.5) = 3.142857$$

$$U\mathbf{x} = \mathbf{z}: \quad x_3 = \frac{3.142857}{1.571429} = 2,\quad x_2 = \frac{11.5 - 0.5(2)}{3.5} = 3,\quad x_1 = \frac{-1 - 3 + 6}{2} = 1$$

> $\mathbf{z}$ es exactamente la columna derecha que quedaba en el Ejemplo 1
> ($-1,\ 23/2,\ 22/7$). Otra vez: es el mismo cálculo.

---

## Ejemplo 5 — LU con pivoteo (el caso del Programa 6-3)

**Columna 1:** el mayor $|a_{i1}|$ es el $3$ (fila 3) → se intercambian filas 1 y 3.
La permutación queda $IP = (3, 2, 1)$.

$$PA = \begin{bmatrix} 3 & 1 & -3 \\ -1 & 3 & 2 \\ 2 & 1 & -3 \end{bmatrix}$$

| Cálculo | Resultado |
|---|---|
| $u_{1j}$ | $3,\ 1,\ -3$ |
| $l_{21} = -1/3,\ \ l_{31} = 2/3$ | $-0.33333,\ 0.66667$ |
| Columna 2: $\lvert 3\rvert > \lvert 1\rvert$ → no se intercambia | — |
| $u_{22} = 3 - (-\tfrac13)(1),\ \ u_{23} = 2 - (-\tfrac13)(-3)$ | $3.33333,\ 1$ |
| $l_{32} = \dfrac{1 - (\tfrac23)(1)}{10/3} = \dfrac{1/3}{10/3}$ | $0.1$ |
| $u_{33} = -3 - (\tfrac23)(-3) - (0.1)(1)$ | $-1.1$ |

$$\det A = (-1)^{1}\cdot(3)(3.3333)(-1.1) = -(-11) = 11 \;\checkmark$$

Un intercambio → signo negativo; el producto de la diagonal da $-11$ y el signo lo corrige a $11$.

---

# PARTE III — PROGRAMA 6-2 PASO A PASO (Inversión)

## Qué hace

> Recibe una matriz cuadrada $A$, arma $[A \mid I]$, aplica eliminación de Gauss con
> pivoteo y **devuelve $A^{-1}$ y $\det A$**.

**Detalle importante:** aunque la teoría habla de Gauss-Jordan, el programa hace
**Gauss + sustitución hacia atrás repetida $N$ veces** (una por cada columna de la
derecha). El resultado es el mismo y se reutiliza la misma rutina del Programa 6-1.

## Estructura del código R

```
machine_eps()    → épsilon de la máquina (el bucle del listado)
gauss_inverse()  → equivalente a SUBROUTINE GAUSS(N,A)
print_mat()      → imprime en el formato 1P5E12.5 del libro
hilbert()        → genera la matriz de Hilbert (problema mal condicionado)
```

## Paso 1 — épsilon de la máquina

```r
eps <- 1.0
while (1.0 + eps > 1.0) eps <- eps / 2
eps * 2
```

Se divide entre 2 hasta que $1 + \varepsilon$ ya no se distingue de $1$. Se usa para
"limpiar" residuos de redondeo (paso 4).

| | Libro (FORTRAN, `REAL*4`) | R (`double`) |
|---|---|---|
| $\varepsilon$ | $5.96\times10^{-8}$ | $2.22\times10^{-16}$ |

> El libro usa precisión simple; R usa doble. Por eso el libro muestra
> `-9.93411E-09` donde debería ir un $0$, y R da el $0$ exacto.

## Paso 2 — arreglo aumentado

```r
M <- cbind(A, diag(N))     # [A | I], N x 2N
```

En el FORTRAN esto está escrito a mano en los `DATA`: `2, 1,-3, 1, 0, 0` etc.

## Paso 3 — pivoteo y determinante

```r
pv <- i
for (j in (i + 1):N)
  if (abs(M[pv, i]) < abs(M[j, i])) pv <- j      # mayor |pivote|
if (pv != i) {
  M[c(i, pv), ] <- M[c(pv, i), ]                 # intercambio (las 2N columnas)
  det <- -det                                    # cada cambio invierte el signo
}
```

Se intercambian las **$2N$ columnas** de la fila, no solo las de $A$. Eso es lo que
hace que la inversa salga bien sin tener que "desordenar" nada al final: intercambiar
dos ecuaciones de $AX = I$ no cambia la solución $X$.

## Paso 4 — eliminación con limpieza de redondeo

```r
r <- M[jr, i] / M[i, i]
for (kc in (i + 1):(2 * N)) {
  temp <- M[jr, kc]
  M[jr, kc] <- M[jr, kc] - r * M[i, kc]
  if (abs(M[jr, kc]) < eps2 * abs(temp)) M[jr, kc] <- 0
}
```

Si el resultado de una resta es menor que $2\varepsilon$ veces el valor anterior, se
considera **ruido de redondeo** y se fuerza a $0$.

> En el FORTRAN la condición es `ABS(A(JR,KC)).LT.EPS2*TEMP` sin `ABS` en `TEMP`, así
> que solo actúa cuando `TEMP` es positivo. En R se usó `abs(temp)` para que funcione
> en ambos signos.

## Paso 5 — determinante

```r
det <- det * prod(diag(M[, 1:N]))
```

Al terminar la eliminación, la parte izquierda es triangular superior: el determinante
es el producto de los pivotes, con el signo que fueron acumulando los intercambios.
**Sale gratis** de la eliminación.

## Paso 6 — sustitución hacia atrás $N$ veces

```r
for (m in (N + 1):(2 * N)) {               # cada columna de la derecha
  M[N, m] <- M[N, m] / M[N, N]
  for (nv in (N - 1):1) {
    va <- M[nv, m] - sum(M[nv, (nv + 1):N] * M[(nv + 1):N, m])
    M[nv, m] <- va / M[nv, nv]
  }
}
```

Es la sustitución hacia atrás de siempre, pero el "lado derecho" es la columna $m$, y
la solución se escribe **encima** de esa misma columna. Al terminar, las columnas
$N+1 \ldots 2N$ son $A^{-1}$.

> **En el FORTRAN:** el bucle `DO 1100 M=N+1,N*2` es exactamente esto.

---

# PARTE IV — PROGRAMA 6-3 PASO A PASO (LU)

## Qué hace

> Factoriza $A = LU$ (o $PA = LU$ con pivoteo), imprime $L$ y $U$ en **forma compacta**,
> el vector de **permutación** y el **determinante**.

## Estructura del código R

```
pivot()      → equivalente a SUBROUTINE PIVOT
lu_decomp()  → el programa principal: renglón de U / columna de L
lu_solve()   → extra: resuelve Ax = y con la factorización (Lz = Py, Ux = z)
show_lu()    → imprime permutación, forma compacta y determinante
```

## La subrutina PIVOT

```r
jj <- j - 1 + which.max(abs(st$A[j:N, j]))     # fila con el mayor |a(k,j)|
if (jj != j) {
  st$A[c(j, jj), ]  <- st$A[c(jj, j), ]        # intercambia en A
  st$EL[c(j, jj), ] <- st$EL[c(jj, j), ]       # y en EL (la parte de L ya hecha)
  st$IP[c(j, jj)]   <- st$IP[c(jj, j)]         # registra la permutación
  st$swaps <- st$swaps + 1
}
```

Tres cosas se intercambian juntas:
1. Las filas de $A$ (para seguir calculando con la matriz permutada).
2. Las filas de `EL`: los $l_{ik}$ ya calculados **viajan con su fila**.
3. El vector `IP`, que recuerda el orden original. Con `IP = (3,2,1)`: la fila 1
   de $PA$ es la fila 3 de $A$, etc.

## El núcleo: renglón de $U$, columna de $L$

```r
for (m in 2:N) {
  if (pivoting) st <- pivot(st, m)

  for (j in m:N) {                                          # renglón m de U
    s <- sum(st$EL[m, 1:(m - 1)] * st$EL[1:(m - 1), j])
    st$EL[m, j] <- st$A[m, j] - s
  }
  for (i in (m + 1):N) {                                    # columna m de L
    s <- sum(st$EL[i, 1:(m - 1)] * st$EL[1:(m - 1), m])
    st$EL[i, m] <- (st$A[i, m] - s) / st$EL[m, m]
  }
}
```

Las dos fórmulas del formulario, literal. Nótese que todo se escribe en la **misma
matriz** `EL`: arriba va $U$, abajo va $L$. Por eso se llama forma compacta.

## El determinante y la rareza de `IPC`

```fortran
      IPC=1                                   ! empieza en 1, no en 0
      ...
      IF (IPC.EQ.INT(IPC/2)*2) DE=-DE         ! si IPC es par, cambia el signo
```

En el FORTRAN, `IPC` **empieza en 1**. Por eso el libro imprime
`NUMERO DE PIVOTEOS = 2` cuando en realidad hubo **un solo** intercambio. La regla
"si `IPC` es par, negar" equivale a "si el número real de intercambios es impar, negar".

En R se usó un contador real (`swaps`, empieza en 0) y la regla directa:

```r
if (st$swaps %% 2 == 1) de <- -de
```

> **Si te preguntan por qué el libro dice 2 pivoteos:** es un detalle de la
> implementación, el contador arranca en 1. Hubo un solo intercambio (filas 1 y 3).

## Detalle sobre el criterio de pivoteo

`PIVOT` compara los $|a_{k,j}|$ de la matriz **original** (permutada), no los valores
ya reducidos $a_{k,j} - \sum l_{k,i}u_{i,j}$. Una LU estrictamente robusta compararía
los reducidos (que son los que realmente terminan en la diagonal). Para las matrices
de este tema ambos criterios eligen la misma fila, así que el resultado coincide.

---

# PARTE V — SALIDAS VERIFICADAS

## Programa 6-2

```
 MATRIZ ORIGINAL
 2.00000E+00 1.00000E+00-3.00000E+00
-1.00000E+00 3.00000E+00 2.00000E+00
 3.00000E+00 1.00000E+00-3.00000E+00

 EPSILON DE LA MAQUINA = 2.2204460e-16

 DETERMINANTE = 11.00000

 MATRIZ INVERSA
-1.00000E+00 0.00000E+00 1.00000E+00
 2.72727E-01 2.72727E-01-9.09091E-02
-9.09091E-01 9.09091E-02 6.36364E-01
```

✅ **Coincide con el libro** (salvo el `-9.93E-09` del libro, que es redondeo de precisión simple).

**Verificaciones extra en R:**

| Prueba | Resultado |
|---|---|
| $\max\lvert AA^{-1} - I\rvert$ | $4.44\times10^{-16}$ |
| Diferencia con `solve(A)` de R | $0$ |
| $A^{-1}\mathbf{y}$ con $\mathbf{y} = (-1, 12, 0)$ | $(1,\ 3,\ 2)$ ✓ |
| Matriz singular (fila 2 = 2·fila 1) | `LA MATRIZ ES SINGULAR` ✓ |

## Programa 6-3 — con pivoteo (el del libro)

```
 Pivoteo en columna 1: se intercambian filas 1 y 3
 PERMUTACION
  3   2   1

 MATRICES LU EN FORMA COMPACTA
 3.00000E+00 1.00000E+00-3.00000E+00
-3.33333E-01 3.33333E+00 1.00000E+00
 6.66667E-01 1.00000E-01-1.10000E+00

 DETERMINANTE = 11.00000
```

✅ **Coincide exactamente con el libro.** Verificación: $\max\lvert LU - PA\rvert = 0$.

## Programa 6-3 — sin pivoteo (el de la teoría)

```
 MATRICES LU EN FORMA COMPACTA
 2.00000E+00 1.00000E+00-3.00000E+00
-5.00000E-01 3.50000E+00 5.00000E-01
 1.50000E+00-1.42857E-01 1.57143E+00

 DETERMINANTE = 11.00000
```

✅ **Coincide con el ejemplo de la teoría.**

> **Dos factorizaciones distintas, mismo determinante.** La LU depende del pivoteo,
> pero $\det A$ no.

## Resolver con LU

| | $\mathbf{z}$ ($L\mathbf{z} = P\mathbf{y}$) | $\mathbf{x}$ ($U\mathbf{x} = \mathbf{z}$) |
|---|---|---|
| Sin pivoteo | $(-1,\ 11.5,\ 3.14286)$ | $(1,\ 3,\ 2)$ |
| Con pivoteo | $(0,\ 12,\ -2.2)$ | $(1,\ 3,\ 2)$ |

Los $\mathbf{z}$ son distintos (porque las $L$ son distintas), pero $\mathbf{x}$ es el mismo.

**Reutilizando la factorización** con $\mathbf{y} = \mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$:

```
   y = (1,0,0)  ->  x = ( -1.00000   0.27273  -0.90909)
   y = (0,1,0)  ->  x = (  0.00000   0.27273   0.09091)
   y = (0,0,1)  ->  x = (  1.00000  -0.09091   0.63636)
```

Son las **columnas de $A^{-1}$**: el Programa 6-3 reproduce el resultado del 6-2.

## El caso que exige pivoteo

Con la matriz del Ejemplo 2 ($a_{11} = 0$):

```
 Sin pivoteo:
 -> PIVOTE CERO (u11 = 0): l21 = a21/0, active el pivoteo

 Con pivoteo:
 PERMUTACION
  3   1   2
 MATRICES LU EN FORMA COMPACTA
 2.00000E+00 4.00000E+00 1.00000E+00
 0.00000E+00 1.00000E+01 1.00000E+00
 5.00000E-01 1.00000E-01-1.60000E+00
 DETERMINANTE = -32.00000

 Con y = (2, 6, 5):  x = (2.7188, 0.4062, -2.0625)
```

Dos intercambios → signo $+$; $2 \cdot 10 \cdot (-1.6) = -32$ ✓

---

# PARTE VI — PROBLEMAS MAL CONDICIONADOS

## La idea

Un sistema **mal condicionado** tiene solución, pero un cambio diminuto en los
coeficientes (o el redondeo) produce un cambio **enorme** en la solución.

Síntomas (del tema):
- $\det(A)\cdot\det(A^{-1})$ se aleja de $1$.
- $(A^{-1})^{-1}$ se aleja de $A$.
- $AA^{-1}$ se aleja de $I$.
- Los elementos de la diagonal son pequeños comparados con los de fuera.

## Geometría en 2D

Dos rectas **casi paralelas**: se cortan, pero un pequeño giro de una de ellas mueve
el punto de corte muchísimo. El caso límite (paralelas exactas) es la matriz singular.

## Ejemplo: matriz de Hilbert $h_{ij} = \dfrac{1}{i+j-1}$

$$H_4 = \begin{bmatrix} 1 & 1/2 & 1/3 & 1/4 \\ 1/2 & 1/3 & 1/4 & 1/5 \\ 1/3 & 1/4 & 1/5 & 1/6 \\ 1/4 & 1/5 & 1/6 & 1/7 \end{bmatrix}
\qquad
H_4^{-1} = \begin{bmatrix} 16 & -120 & 240 & -140 \\ -120 & 1200 & -2700 & 1680 \\ 240 & -2700 & 6480 & -4200 \\ -140 & 1680 & -4200 & 2800 \end{bmatrix}$$

Elementos de $H$ menores que $1$; elementos de $H^{-1}$ de hasta $6480$. Ese contraste
ya es una señal.

**Ejercicio propuesto, resuelto en R:**

| Cálculo | Ideal | Obtenido |
|---|:---:|---|
| $\det H_4$ | — | $1.653\times10^{-7}$ |
| $\det H_4^{-1}$ | — | $6.048\times10^{6}$ |
| b) $\det H \cdot \det H^{-1}$ | $1$ | $0.999999999999807$ |
| a) $\max\lvert H^{-1}(H^{-1})^{-1} - I\rvert$ | $0$ | $2.27\times10^{-13}$ |

En doble precisión, $H_4$ todavía se defiende (pierde ~3 cifras de 16). El problema
aparece al crecer $n$:

| $n$ | $\kappa(H)$ | $\max\lvert HH^{-1} - I\rvert$ | $\det\cdot\det(\text{inv}) - 1$ |
|:---:|---:|---:|---:|
| 3 | $5.2\times10^{2}$ | $7.1\times10^{-15}$ | $-3.1\times10^{-15}$ |
| 4 | $1.6\times10^{4}$ | $1.1\times10^{-13}$ | $-1.9\times10^{-13}$ |
| 6 | $1.5\times10^{7}$ | $1.5\times10^{-10}$ | $1.3\times10^{-10}$ |
| 8 | $1.5\times10^{10}$ | $6.0\times10^{-8}$ | $1.8\times10^{-9}$ |
| 10 | $1.6\times10^{13}$ | $8.0\times10^{-5}$ | $-2.2\times10^{-5}$ |
| 12 | $1.8\times10^{16}$ | **$8.8\times10^{-2}$** | **$-8.9\times10^{-3}$** |

$\kappa(A) = \lVert A\rVert\,\lVert A^{-1}\rVert$ es el **número de condición**.

> **Regla práctica:** se pierden aproximadamente $\log_{10}\kappa(A)$ cifras
> significativas. Con $\kappa \approx 10^{16}$ y 16 cifras disponibles, no queda
> **ninguna**: con $n = 12$, $HH^{-1}$ difiere de $I$ en casi un 9 %.

> En el libro (precisión simple, ~7 cifras) $H_4$ con $\kappa \approx 1.6\times10^4$ ya
> pierde más de la mitad de las cifras. Por eso allí el ejemplo se ve mucho peor.

---

# PARTE VII — GUION DE LA EXPOSICIÓN

## Lámina 1 — Portada
> "Capítulo 6 de Nakamura: Álgebra Lineal Numérica. La teoría ya la vimos en Cálculo;
> lo que vamos a ver hoy es cómo se lleva a código y qué problemas aparecen cuando la
> computadora trabaja con un número finito de cifras. Vamos a mostrar dos programas: el
> 6-2, que invierte una matriz, y el 6-3, que hace la descomposición LU."

## Lámina 2 — Motivación
> "Resolver $A\mathbf{x} = \mathbf{y}$ está detrás de casi todo lo que vimos: los
> splines del tema 2, el sistema de Vandermonde del tema 5, mínimos cuadrados,
> ecuaciones diferenciales. ¿Y por qué no usar Cramer? Porque cuesta del orden de $N!$.
> Con 20 ecuaciones son $10^{18}$ operaciones; Gauss lo hace en unas 5 000."

## Lámina 3 — Eliminación de Gauss
> "Dos fases: eliminar hacia adelante hasta dejar una triangular superior, y despejar
> de abajo hacia arriba. Con este sistema sale $x = (1, 3, 2)$. Guarden esta matriz,
> porque la vamos a usar en todos los ejemplos."

## Lámina 4 — Pivoteo
> "Gauss divide por el pivote. Si el pivote es cero, el método se cae; si es muy
> pequeño, amplifica el error de redondeo. El pivoteo parcial busca en la columna el
> elemento de mayor valor absoluto y lo sube. Intercambiar ecuaciones no cambia la
> solución, pero cada intercambio cambia el signo del determinante."

## Lámina 5 — Inversa por Gauss-Jordan
> "Para invertir, se pega la identidad a la derecha y se hacen las mismas operaciones
> sobre todo. Cuando a la izquierda queda $I$, a la derecha queda $A^{-1}$. La razón:
> cada operación de fila es multiplicar por una matriz; si el producto de todas
> convierte $A$ en $I$, ese producto es $A^{-1}$, y eso es lo que le pasó a la
> identidad de la derecha."

## Lámina 6 — Programa 6-2: estructura
> "El programa arma $[A|I]$, que es una matriz de $N \times 2N$, y hace eliminación de
> Gauss sobre las $2N$ columnas. Después aplica la sustitución hacia atrás $N$ veces,
> una por cada columna de la derecha. Técnicamente es Gauss y no Gauss-Jordan, pero
> el resultado es el mismo."

## Lámina 7 — Programa 6-2: detalles
> "Tres detalles del listado. Primero, calcula el épsilon de la máquina dividiendo por
> 2 hasta que $1 + \varepsilon$ ya no se distingue de 1. Segundo, usa ese épsilon para
> poner en cero los residuos de redondeo. Tercero, calcula el determinante de paso:
> producto de los pivotes, con un cambio de signo por cada intercambio."

## Lámina 8 — Programa 6-2: salida
> "Reproduce la salida del libro. La única diferencia es este $-9.9\times10^{-9}$ que
> en el libro aparece donde debería haber un cero: es redondeo de precisión simple. En
> R, con doble precisión, sale el cero exacto. Si multiplicamos la inversa por 11, que
> es el determinante, aparecen enteros: es la adjunta. Y verificamos que
> $AA^{-1} = I$ con error de $10^{-16}$."

## Lámina 9 — Descomposición LU
> "La idea: escribir $A$ como producto de una triangular inferior con unos en la
> diagonal y una triangular superior. Se calcula alternando: un renglón de $U$, una
> columna de $L$, otro renglón de $U$... Cada elemento usa solo los que ya se
> calcularon."

## Lámina 10 — Gauss y LU son lo mismo
> "Lo interesante es que $L$ son exactamente los multiplicadores que usamos en Gauss:
> $-1/2$, $3/2$ y $-1/7$. Y $U$ es la triangular que quedó al final. LU es Gauss
> guardando los multiplicadores. Como sobra espacio, se guardan juntas en una sola
> matriz: la forma compacta."

## Lámina 11 — Programa 6-3
> "El programa sigue las fórmulas al pie de la letra. Con pivoteo, la subrutina PIVOT
> intercambia tres cosas a la vez: las filas de $A$, las filas de $L$ que ya se
> calcularon, y un vector IP que recuerda la permutación. Así lo que se factoriza es
> $PA$."

## Lámina 12 — Programa 6-3: salida
> "Con pivoteo sale la permutación 3, 2, 1 y la forma compacta del libro. Un detalle:
> el libro dice 'número de pivoteos = 2', pero hubo uno solo. El contador del FORTRAN
> arranca en 1. Sin pivoteo sale la LU de la teoría. Dos factorizaciones distintas,
> y el mismo determinante: 11."

## Lámina 13 — Resolver con LU
> "¿Para qué sirve? Para resolver $A\mathbf{x} = \mathbf{y}$ en dos pasos triangulares:
> $L\mathbf{z} = P\mathbf{y}$ hacia adelante y $U\mathbf{x} = \mathbf{z}$ hacia atrás.
> La ventaja es que factorizar cuesta $N^3$, pero cada nuevo $\mathbf{y}$ cuesta solo
> $N^2$. Si le paso las columnas de la identidad, obtengo las columnas de $A^{-1}$: el
> mismo resultado del Programa 6-2."

## Lámina 14 — Determinantes
> "El determinante por definición tiene $N!$ términos. Pero como el determinante de un
> producto es el producto de determinantes, y el de una triangular es el producto de
> su diagonal, $\det A = \det U$. En nuestro caso, $2 \cdot 3.5 \cdot 1.5714 = 11$.
> Sale gratis de la factorización."

## Lámina 15 — Mal condicionados
> "Y la advertencia del capítulo. Hay matrices que tienen inversa, pero el redondeo la
> destruye. El ejemplo clásico es la de Hilbert. Medimos $HH^{-1} - I$, que debería
> ser cero: con $n = 4$ el error es $10^{-13}$, pero con $n = 12$ es casi un 9 %. La
> regla es que se pierden tantas cifras como el logaritmo del número de condición."

## Lámina 16 — Conclusiones
> "Resumiendo: Gauss es la base de todo; el pivoteo no es opcional; la inversa sale de
> resolver $N$ sistemas con la identidad; LU es Gauss guardando los multiplicadores y
> permite reutilizar la factorización; el determinante sale gratis; y hay que vigilar
> el condicionamiento porque el álgebra exacta y la aritmética de la computadora no
> siempre coinciden."

---

# PARTE VIII — PREGUNTAS PROBABLES EN LA DEFENSA

**¿Por qué no se resuelve $A\mathbf{x} = \mathbf{y}$ calculando $A^{-1}\mathbf{y}$?**
Porque invertir cuesta ~$2N^3$ y además acumula más redondeo. Con LU se gasta
$\tfrac23N^3$ una vez y $2N^2$ por cada $\mathbf{y}$. La inversa solo vale la pena
cuando se necesita la matriz en sí.

**¿Qué pasa si no se pivotea?**
Si un pivote es $0$, división entre cero. Si es muy pequeño, el multiplicador $r$ es
enorme y amplifica los errores de redondeo. El pivoteo parcial hace que
$\lvert r\rvert \le 1$ siempre.

**¿Por qué el pivoteo cambia el signo del determinante?**
Intercambiar dos filas de una matriz cambia el signo del determinante. Con $s$
intercambios, el factor es $(-1)^s$.

**¿La LU es única?**
Con la convención de unos en la diagonal de $L$ (Doolittle) y un orden de filas
fijo, sí. Pero con pivoteo y sin pivoteo salen LU distintas, porque se factorizan
matrices distintas ($A$ y $PA$). El determinante sí es el mismo.

**¿Cómo detecta el programa una matriz singular?**
Si después de pivotear el mejor pivote disponible es $0$, toda la columna es cero
por debajo de la diagonal y no hay forma de continuar: `LA MATRIZ ES SINGULAR`.

**¿Para qué sirven las dos matrices `A` y `EL` en el 6-3?**
Solo para la demostración. El libro dice que se puede eliminar una y escribir todo en
`A`, porque cada $a_{ij}$ se usa una sola vez, justo antes de reemplazarlo por $l_{ij}$
o $u_{ij}$.

**¿Por qué el libro dice "2 pivoteos" si hubo uno?**
El contador `IPC` arranca en 1. La regla del signo está ajustada a eso.

**¿Qué es un problema mal condicionado y cómo se reconoce?**
Uno donde un cambio pequeño en los datos produce un cambio grande en la solución. Se
reconoce con $\kappa(A)$ grande, o en la práctica comprobando $AA^{-1} \approx I$ y
$\det A\cdot\det A^{-1} \approx 1$.

---

# PARTE IX — RESUMEN COMPARATIVO

| | Gauss | Gauss-Jordan / Inversa | LU |
|---|---|---|---|
| Resultado | $\mathbf{x}$ | $A^{-1}$ | $L$, $U$, $P$ |
| Costo | $\tfrac23N^3$ | $\approx 2N^3$ | $\tfrac23N^3$ + $2N^2$ por $\mathbf{y}$ |
| Varios $\mathbf{y}$ | repetir todo | multiplicar $A^{-1}\mathbf{y}$ | solo sustituciones |
| Determinante | gratis | gratis | gratis |
| Programa | 6-1 | **6-2** | **6-3** |

## Decisión rápida

```
¿Qué necesito?
├── Resolver un solo sistema         → GAUSS con pivoteo
├── Resolver muchos con la misma A   → LU (factorizar una vez)
├── La matriz inversa en sí          → [A | I] (Programa 6-2)
└── El determinante                  → producto de pivotes (nunca la definición)

¿Algún pivote es 0 o muy pequeño?
└── Siempre                          → PIVOTEO PARCIAL

¿Los resultados no cuadran (A·A⁻¹ ≠ I)?
└── Revisar el condicionamiento κ(A)
```

> **Frase para cerrar:** *"LU es Gauss con memoria: guarda los multiplicadores para no
> tener que repetir la eliminación cada vez que cambia el lado derecho."*
