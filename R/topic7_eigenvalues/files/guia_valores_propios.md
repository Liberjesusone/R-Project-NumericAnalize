# Guía de Valores Propios — Capítulo 7 (Nakamura)

> Universidad de Los Andes · CeSiMo · Análisis Numérico
> Programas 7-1 (Método de interpolación) y 7-2 (Householder/bisección) implementados en R

---

# PARTE I — FORMULARIO

## 1. Definición

$$A\mathbf{x} = \lambda\mathbf{x} \quad\Longleftrightarrow\quad (A - \lambda I)\mathbf{x} = \mathbf{0}$$

Tiene solución no trivial ($\mathbf{x}\ne\mathbf{0}$) **solo si**

$$\boxed{\;f(\lambda) \equiv \det(A - \lambda I) = 0\;}$$

$f(\lambda)$ es el **polinomio característico**, de grado $N$ → hay $N$ valores propios
(contando multiplicidad, posiblemente complejos).

**Chequeos rápidos:**

$$\sum_i \lambda_i = \operatorname{traza}(A), \qquad \prod_i \lambda_i = \det A = f(0)$$

**Propiedades:**
- $A$ simétrica → todos los $\lambda$ son **reales** (y si todos son positivos, $A$ es *definida positiva*).
- $A$ real no simétrica → puede tener $\lambda$ complejos, siempre en **pares conjugados**.
- Una **transformación de similaridad** $B = PAP^{-1}$ **no cambia** los valores propios.

---

## 2. Método de interpolación (Programa 7-1)

1. Evaluar $f(\lambda_j) = \det(A - \lambda_j I)$ en $\lambda_j = j\,\Delta\lambda$, $j = 0,\dots,N$.
2. Tabla de diferencias hacia adelante: $\Delta^k f_j = \Delta^{k-1}f_{j+1} - \Delta^{k-1}f_j$.
3. Fórmula de Newton hacia adelante, con $s = \lambda/\Delta\lambda$:

$$f = \sum_{k=0}^{N} \Delta^k f_0 \binom{s}{k}, \qquad \binom{s}{k} = \frac{s(s-1)\cdots(s-k+1)}{k!}$$

4. **Coeficientes de Markov** = coeficientes de $\binom{s}{k}$ en potencias de $s$.
5. Cambio de variable: $c_i^{(\lambda)} = c_i^{(s)} / \Delta\lambda^{\,i}$.

> Como $f$ es un polinomio de grado $N$, **$N+1$ puntos lo determinan exactamente**.
> No hay error de interpolación, solo de redondeo.

---

## 3. Transformación de Householder

$$A^{(m+1)} = P\,A^{(m)}\,P, \qquad P = I - \frac{\mathbf{u}\mathbf{u}^T}{h}$$

$$G = \operatorname{signo}(a_{m+1,m})\left[\sum_{i=m+1}^{N} a_{i,m}^2\right]^{1/2}, \qquad h = G^2 + G\,a_{m+1,m}$$

$$\mathbf{u} = (0,\dots,0,\; a_{m+1,m} + G,\; a_{m+2,m},\; \dots,\; a_{N,m})^T$$

**Propiedades de $P$:** $P^T = P$, $\;PP = I$, $\;P^{-1} = P$ → es una similaridad.

Con $N-2$ pasos, una matriz **simétrica** queda **tridiagonal**.

---

## 4. Determinante de una tridiagonal (ec. 7.3.15)

$$M = \begin{bmatrix} a_1 & b_1 & & \\ b_1 & a_2 & b_2 & \\ & b_2 & a_3 & \ddots \\ & & \ddots & a_N \end{bmatrix}$$

$$p_0 = 1, \qquad p_1 = a_1 - \lambda, \qquad \boxed{\;p_i = (a_i - \lambda)\,p_{i-1} - b_{i-1}^2\,p_{i-2}\;}$$

$p_i(\lambda)$ es el polinomio característico de la **submatriz principal** $i\times i$.
Cuesta $O(N)$ en lugar de $O(N^3)$.

**Entrelazamiento:** los valores propios de la submatriz $(L-1)$ **separan** a los de la submatriz $L$:

$$\lambda_1^{(L)} < \lambda_1^{(L-1)} < \lambda_2^{(L)} < \lambda_2^{(L-1)} < \cdots < \lambda_{L-1}^{(L-1)} < \lambda_L^{(L)}$$

→ cada intervalo contiene **exactamente una** raíz de $p_L$.

---

## 5. Método de las potencias

$$\mathbf{v} = A\mathbf{u}^{(k)}, \qquad \lambda^{(k+1)} = \frac{(\mathbf{v},\mathbf{v})}{(\mathbf{v},\mathbf{u}^{(k)})}, \qquad \mathbf{u}^{(k+1)} = \frac{\mathbf{v}}{\lambda^{(k+1)}}$$

Converge al valor propio de **mayor $|\lambda|$**, con error $\sim |\lambda_{N-1}/\lambda_N|^k$.

**Potencias inversas:** igual, pero resolviendo $A\mathbf{v} = \mathbf{u}^{(k)}$ (es aplicar
potencias a $A^{-1}$) → converge al de **menor $|\lambda|$**.

---

## 6. Iteración QR

$$A_k = Q_k R_k \;\Rightarrow\; A_{k+1} = R_k Q_k = Q_k^{-1} A_k Q_k$$

Con desplazamiento: $A_k - \alpha_k I = Q_kR_k$, $\;A_{k+1} = R_kQ_k + \alpha_k I$.

Converge a una **triangular superior por bloques**: bloques $1\times1$ (valores reales) y
$2\times2$ (pares complejos).

**Raíces de un polinomio** $a_N x^N + \cdots + a_0$ → valores propios de la **matriz de Frobenius**:

$$P = \begin{bmatrix} -p_{N-1} & -p_{N-2} & \cdots & -p_0 \\ 1 & 0 & \cdots & 0 \\ & \ddots & & \vdots \\ 0 & & 1 & 0 \end{bmatrix}, \qquad p_i = \frac{a_i}{a_N}$$

---

# PARTE II — TEORÍA CON EJEMPLOS RESUELTOS

## Ejemplo 1 — El sistema homogéneo del tema

$$3x - 2y + z = 0, \qquad x + y + 2z = 0, \qquad 4x - y + 3z = 0$$

$\det A = 0$ (la tercera ecuación es la suma de las dos primeras), así que hay soluciones
no triviales: $x = \beta,\; y = \beta,\; z = -\beta$.

> **La conexión:** pedir $A\mathbf{x} = \lambda\mathbf{x}$ con $\mathbf{x}\ne\mathbf{0}$ es pedir
> que el sistema homogéneo $(A-\lambda I)\mathbf{x} = \mathbf{0}$ tenga solución no trivial.
> Por eso se busca el $\lambda$ que hace $\det(A - \lambda I) = 0$.

---

## Ejemplo 2 — Aplicación física: masas y resortes

Tres masas iguales unidas por resortes. Proponiendo $y_j = e^{i\omega t} f_j$, la ecuación
$m y'' = -K y$ se convierte en

$$\frac{K}{m}\mathbf{f} = \omega^2\mathbf{f} \quad\Rightarrow\quad \lambda = \omega^2$$

Con todos los $k = 1$ y $m = 1$:

$$A = \begin{bmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{bmatrix}, \qquad
\lambda = 2-\sqrt2,\; 2,\; 2+\sqrt2, \qquad \omega = 0.765,\; 1.414,\; 1.848 \text{ rad/s}$$

| Modo | $\omega$ | Vector propio | Movimiento |
|:---:|:---:|---|---|
| 1 | 0.765 | $(0.5,\ 0.707,\ 0.5)$ | las tres masas en fase |
| 2 | 1.414 | $(-0.707,\ 0,\ 0.707)$ | la del centro quieta |
| 3 | 1.848 | $(-0.5,\ 0.707,\ -0.5)$ | alternadas |

> **Los valores propios son las frecuencias naturales** y los vectores propios son las
> formas de vibrar. Es la matriz simétrica y tridiagonal típica de la ingeniería.

> En la transcripción del tema la matriz aparece con los elementos fuera de la diagonal
> positivos; físicamente van con signo **negativo** ($-k/m$), como arriba.

---

## Ejemplo 3 — Polinomio característico por interpolación ($N = 4$, $\Delta\lambda = 0.5$)

$$A = \begin{bmatrix} 4 & 3 & 2 & 1 \\ 3 & 3 & 2 & 1 \\ 2 & 2 & 2 & 1 \\ 1 & 1 & 1 & 1 \end{bmatrix}$$

**Paso 1 — 5 determinantes:**

| $\lambda$ | 0 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|:---:|:---:|:---:|:---:|:---:|
| $f(\lambda)$ | 1 | 0.0625 | 0 | −4.4375 | −17 |

**Paso 2 — tabla de diferencias** (primera fila: $\Delta^k f_0$):

| $\lambda$ | $f$ | $\Delta$ | $\Delta^2$ | $\Delta^3$ | $\Delta^4$ |
|---|---:|---:|---:|---:|---:|
| 0.0 | 1 | **−0.9375** | **0.875** | **−5.25** | **1.5** |
| 0.5 | 0.0625 | −0.0625 | −4.375 | −3.75 | |
| 1.0 | 0 | −4.4375 | −8.125 | | |
| 1.5 | −4.4375 | −12.5625 | | | |
| 2.0 | −17 | | | | |

**Paso 3 — Newton hacia adelante** ($s = \lambda/0.5 = 2\lambda$):

$$f = 1 - 0.9375\,s + 0.875\binom{s}{2} - 5.25\binom{s}{3} + 1.5\binom{s}{4}$$

**Paso 4 — coeficientes de Markov** (los binomiales en potencias de $s$):

| $k$ | $\binom{s}{k}$ desarrollado | Markov |
|:---:|---|---|
| 2 | $\frac{s^2 - s}{2}$ | $-0.5,\ 0.5$ |
| 3 | $\frac{s^3 - 3s^2 + 2s}{6}$ | $0.3333,\ -0.5,\ 0.1667$ |
| 4 | $\frac{s^4 - 6s^3 + 11s^2 - 6s}{24}$ | $-0.25,\ 0.4583,\ -0.25,\ 0.0417$ |

Sumando potencia por potencia:

$$f = 1 - 3.5\,s + 3.75\,s^2 - 1.25\,s^3 + 0.0625\,s^4$$

(por ejemplo, el de $s$: $-0.9375 - 0.4375 - 1.75 - 0.375 = -3.5$)

**Paso 5 — de $s$ a $\lambda$** (dividir $c_i$ entre $0.5^i$):

$$\boxed{\;f(\lambda) = \lambda^4 - 10\lambda^3 + 15\lambda^2 - 7\lambda + 1\;}$$

**Chequeos:** $c_0 = 1 = \det A$ ✓, $\;-c_3 = 10 = $ traza ✓, $\;f(1) = 1 - 10 + 15 - 7 + 1 = 0$ → $\lambda = 1$ es valor propio ✓.

> **Detalle bonito:** $\Delta^4 f = 1.5 = 4!\,(0.5)^4$. La $N$-ésima diferencia de un
> polinomio de grado $N$ es constante: $N!\,h^N$ por el coeficiente principal. Es lo
> mismo que vimos en diferenciación: $\Delta^N/h^N$ aproxima la derivada $N$-ésima.

---

## Ejemplo 4 — Householder a mano (ejemplo del tema)

$$A = \begin{bmatrix} 3 & 1 & 4 \\ 1 & 7 & 2 \\ 4 & 2 & 0 \end{bmatrix}$$

Solo hace falta **un** paso ($N - 2 = 1$), sobre la columna 1:

$$G = +\sqrt{1^2 + 4^2} = \sqrt{17} = 4.1231, \qquad h = 17 + 4.1231\cdot1 = 21.1231$$

$$\mathbf{u} = (0,\; 1 + 4.1231,\; 4) = (0,\; 5.1231,\; 4)$$

$$P = I - \frac{\mathbf{u}\mathbf{u}^T}{h} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & -0.2425 & -0.9701 \\ 0 & -0.9701 & 0.2425 \end{bmatrix}
= \begin{bmatrix} 1 & 0 & 0 \\ 0 & -1/\sqrt{17} & -4/\sqrt{17} \\ 0 & -4/\sqrt{17} & 1/\sqrt{17} \end{bmatrix}$$

**Fila 1 de $PAP$:** $\;\bigl(3,\ \ 1(-0.2425) + 4(-0.9701),\ \ 1(-0.9701) + 4(0.2425)\bigr) = (3,\ -4.1231,\ 0)$ ✓

$$A^{(2)} = PAP = \begin{bmatrix} 3 & -4.1231 & 0 \\ -4.1231 & 1.3529 & 3.4118 \\ 0 & 3.4118 & 5.6471 \end{bmatrix}$$

**¿Por qué $PP = I$?** Porque $\mathbf{u}^T\mathbf{u} = 2h$:

$$\mathbf{u}^T\mathbf{u} = (a+G)^2 + \sum_{i>m+1} a_{i}^2 = a^2 + 2aG + G^2 + (G^2 - a^2) = 2(G^2 + aG) = 2h$$

y entonces $PP = I - 2\frac{\mathbf{u}\mathbf{u}^T}{h} + \frac{\mathbf{u}(\mathbf{u}^T\mathbf{u})\mathbf{u}^T}{h^2} = I$.

> **¿Por qué el signo de $G$ igual al de $a_{m+1,m}$?** Para que $a + G$ sea una **suma**
> de números del mismo signo y no una resta de números parecidos (cancelación).
> Es la misma idea de siempre: evitar restar números casi iguales.

---

## Ejemplo 5 — Valores propios por la recurrencia

Con la tridiagonal anterior ($a = 3,\ 1.3529,\ 5.6471$; $b = -4.1231,\ 3.4118$):

| Submatriz | Polinomio | Raíces |
|:---:|---|---|
| $1\times1$ | $p_1 = 3 - \lambda$ | $3$ |
| $2\times2$ | $p_2 = (1.3529 - \lambda)p_1 - 17$ | $-2.0281,\ 6.3810$ |
| $3\times3$ | $p_3 = (5.6471 - \lambda)p_2 - 11.640\,p_1$ | $-2.8941,\ 4.3861,\ 8.5080$ |

Entrelazamiento: $-2.0281$ cae entre $-2.8941$ y $4.3861$; $6.3810$ entre $4.3861$ y $8.5080$;
y $3$ entre $-2.0281$ y $6.3810$. ✓

> El tema dice $-2.0289$; el valor correcto es $-2.0281$ (comprobación: la suma
> $-2.0281 + 6.3810 = 4.3529$ debe ser la traza $3 + 1.3529$). Es un error de transcripción.

---

## Ejemplo 6 — Método de las potencias (tabla 7.3)

Misma $A$ 3×3, $\mathbf{u}^{(0)} = (1,1,1)$:

$$\mathbf{v} = A\mathbf{u}^{(0)} = (8,\ 10,\ 6), \qquad \lambda^{(1)} = \frac{64+100+36}{8+10+6} = \frac{200}{24} = 8.3333$$

$$\mathbf{u}^{(1)} = \frac{(8, 10, 6)}{8.3333} = (0.96,\ 1.20,\ 0.72)$$

| $k$ | $\lambda^{(k)}$ | $\mathbf{u}^{(k)}$ |
|:---:|---|---|
| 1 | 8.333333 | 0.9600000, 1.2000000, 0.7200000 |
| 2 | 8.453461 | 0.8233315, 1.2775833, 0.7381592 |
| 5 | 8.506876 | 0.7462269, 1.3558519, 0.6724730 |
| 10 | 8.507999 | 0.7319253, 1.3669517, 0.6655880 |
| 20 | 8.508000 | 0.7314031, 1.3673721, 0.6652983 |

> En la tabla transcrita del tema, el orden de los componentes de $\mathbf{u}$ está mezclado
> en varias filas. El componente **mayor** (~1.367) es siempre el **segundo**.

Converge a $\lambda_{\max} = 8.508$. La rapidez la da $|\lambda_2/\lambda_1| = 4.386/8.508 = 0.516$.

---

# PARTE III — PROGRAMA 7-1 PASO A PASO (Interpolación)

## Qué hace (y qué no hace)

> **No da los valores propios directamente.** Recibe una matriz y un $\Delta\lambda$ y
> devuelve los **coeficientes del polinomio característico**. Las raíces se buscan después
> (en R lo hicimos con Frobenius + QR, la sección 4 del tema).

**Ventaja:** sirve para **cualquier** matriz, simétrica o no, sin desarrollar el
determinante simbólicamente.

## Estructura del código R

```
detmnt()           → SUBROUTINE DETMNT: determinante por Gauss con pivoteo (cap. 6)
charpoly_interp()  → programa principal: N+1 dets, tabla de diferencias, Markov
poly_from_roots()  → coeficientes exactos desde eigen(), para verificar
frobenius()        → matriz de Frobenius de un polinomio
qr_eigen()         → iteración QR sin desplazamiento, lee bloques 1x1 y 2x2
```

## Paso 1 — los $N+1$ determinantes

```r
RA <- (0:N) * DE
for (jj in 0:N)
  FF[jj + 1] <- detmnt(B - RA[jj + 1] * diag(N))
```

`B` se guarda intacta y en cada paso se resta $\lambda$ de la diagonal. El FORTRAN copia
`B` en `A` y hace `A(J,J)=A(J,J)-RA(JJ)` porque `DETMNT` destruye la matriz que recibe.

`detmnt()` es la eliminación de Gauss del capítulo 6: producto de pivotes y cambio de
signo por cada intercambio.

## Paso 2 — tabla de diferencias

```r
DF[, 1] <- FF
for (j in 1:N)
  for (i in 0:(N - j))
    DF[i + 1, j + 1] <- DF[i + 2, j] - DF[i + 1, j]
```

Es la tabla triangular del tema 2 (interpolación de Newton). Solo se usa la **primera
fila** ($\Delta^k f_0$).

## Paso 3 — coeficientes de Markov

```r
MV[2] <- 1                        # C(s,1) = s
CC[1] <- DF[1, 1]                 # Delta^0 f_0
CC[2] <- DF[1, 2]                 # Delta^1 f_0 * s
for (K in 2:N)
  for (L in K:1) {
    MV[L + 1] <- (MV[L] - (K - 1) * MV[L + 1]) / K
    CC[L + 1] <- CC[L + 1] + MV[L + 1] * DF[1, K + 1]
  }
```

**Qué hace la línea de `MV`:** pasa de $\binom{s}{K-1}$ a $\binom{s}{K}$ multiplicando por
$\frac{s - (K-1)}{K}$. Multiplicar un polinomio por $s$ **corre** los coeficientes una
posición (`MV[L]` → posición $L+1$); multiplicar por $-(K-1)$ los **escala**. Por eso:

$$\text{nuevo}_L = \frac{\text{viejo}_{L-1} - (K-1)\,\text{viejo}_L}{K}$$

El bucle va de $L = K$ **hacia abajo** para no pisar un valor que todavía se necesita
(como en la sustitución hacia atrás).

**La línea de `CC`** acumula $\Delta^K f_0 \cdot \binom{s}{K}$ potencia por potencia.

## Paso 4 — cambio de variable

```r
CC <- CC / DE^(0:N)
```

Como $s = \lambda/\Delta\lambda$, se tiene $s^i = \lambda^i/\Delta\lambda^i$.

## Rarezas del listado transcrito

- El `PRINT 223` de los coeficientes de Markov y un `END DO` quedaron fuera de lugar.
  Según la salida del libro, se imprimen dentro del bucle de `K` (una fila por cada $K$),
  y así se hizo en R.
- `DETMNT` imprime `LA MATRIZ ES SINGULAR` y **no asigna** `DET` si $a_{NN} = 0$. Justo
  cuando $\lambda$ es valor propio, el determinante debería ser $0$. En R se devuelve $0$.

---

# PARTE IV — PROGRAMA 7-2 PASO A PASO (Householder/bisección)

## Qué hace

> Solo para matrices **simétricas**. Etapa 1: reduce $A$ a tridiagonal con Householder.
> Etapa 2: calcula los valores propios de **todas las submatrices principales**, de
> $1\times1$ hasta $N\times N$, con bisección + interpolación lineal, usando el
> entrelazamiento para encerrar cada raíz.

## Estructura del código R

```
householder_tridiag() → etapa 1 (bucle DO IR=1,N-2 del programa principal)
tri_det()             → SUBROUTINE DETERM: recurrencia p_i
bisec()               → SUBROUTINE BISEC: bisección + falsa posición
tridiag_eigen()       → etapa 2: recorre L = 1..N con los intervalos entrelazados
power_method()        → extra: potencias directas e inversas (sección 3 del tema)
```

## Etapa 1 — Householder

```r
s   <- sum(A[idx, ir]^2)                # G^2
w   <- if (A[ir + 1, ir] < 0) -1 else 1 # signo(a_{ir+1,ir})
ssr <- sqrt(s)
h   <- s + abs(A[ir + 1, ir]) * ssr     # h = G^2 + G a
u[idx]    <- A[idx, ir]
u[ir + 1] <- A[ir + 1, ir] + ssr * w    # a + G

Au  <- as.vector(A %*% u)
uAu <- sum(u * Au)
A <- A - (outer(u, Au) + outer(Au, u)) / h + (uAu / h^2) * outer(u, u)
```

**La última línea es $PAP$ sin construir $P$.** Desarrollando:

$$PAP = \Bigl(I - \tfrac{\mathbf{u}\mathbf{u}^T}{h}\Bigr)A\Bigl(I - \tfrac{\mathbf{u}\mathbf{u}^T}{h}\Bigr)
= A - \frac{\mathbf{u}(A\mathbf{u})^T + (A\mathbf{u})\mathbf{u}^T}{h} + \frac{\mathbf{u}^TA\mathbf{u}}{h^2}\,\mathbf{u}\mathbf{u}^T$$

> **En el FORTRAN:** `B23` acumula $-(u_i(A\mathbf{u})_j + (A\mathbf{u})_iu_j)$, la matriz
> `T` guarda $A + B23/H$, `UAU` es $\mathbf{u}^TA\mathbf{u}/h^2$, y al final
> `A(I,J)=T(I,J)+UAU*U(I)*U(J)`. Es exactamente la fórmula de arriba.

Después se ponen en cero los elementos con $|a_{ij}| < 10^{-6}$ (umbral del libro). Es
absoluto: con una matriz de elementos muy pequeños borraría información real.

## Etapa 2 — la recurrencia (`DETERM`)

```r
p0 <- 1
p1 <- a[1] - x
for (k in 2:L) {
  p2 <- (a[k] - x) * p1 - b[k - 1]^2 * p0
  p0 <- p1; p1 <- p2
}
```

Evalúa $\det(T_L - xI)$ de la submatriz $L\times L$ en $O(L)$ operaciones. Como la
bisección necesita **muchas** evaluaciones, esta es la razón de tridiagonalizar primero.

## Etapa 2 — los intervalos

```r
ext <- c(-R, EI[[L - 1]], R)                 # los de L-1 separan a los de L
EI[[L]] <- sapply(1:L, function(j) bisec(a, b, L, ext[j], ext[j + 1]))
```

Para la submatriz $L$, la raíz $j$ está entre la raíz $j-1$ y la raíz $j$ de la submatriz
$L-1$. Los extremos exteriores son $\pm R$.

> **En el FORTRAN** los extremos son `EI(L,0)=-99` y `EI(L,L+1)=99`: solo funciona si
> todos los valores propios están en $[-99, 99]$. En R se usó la **cota de Gershgorin**
> $R = \max_i\bigl(|a_i| + |b_{i-1}| + |b_i|\bigr)$, que vale para cualquier matriz.

## Etapa 2 — `BISEC`

```r
if (xh - xl > 1) {
  xm <- (xl + xh) / 2                        # bisección: segura pero lenta
} else {
  xm <- (xl * yh - xh * yl) / (yh - yl)      # interpolación lineal: más rápida
}
```

Mientras el intervalo es grande, bisección. Cuando ya es menor que 1, falsa posición
(la recta entre $(x_l, y_l)$ y $(x_h, y_h)$ corta el eje en $x_m$). Es el tema 3
(ecuaciones no lineales) aplicado a $p_L(\lambda) = 0$.

> **El listado transcrito está incompleto:** la rama de bisección no actualiza `XL`/`XH`
> y hay un `GOTO 80` sin la etiqueta 80. En R se implementó el algoritmo que describe el
> texto, con una protección: si el redondeo deja $x_m$ fuera del intervalo, se usa el
> punto medio.

---

# PARTE V — SALIDAS VERIFICADAS

## Programa 7-1

```
 TABLA DE DIFERENCIAS PARA LOS DETERMINANTES
  0.0000   1.000E+00  -9.375E-01   8.750E-01  -5.250E+00   1.500E+00
  0.5000   6.250E-02  -6.250E-02  -4.375E+00  -3.750E+00
  1.0000   1.665E-16  -4.437E+00  -8.125E+00
  1.5000  -4.437E+00  -1.256E+01
  2.0000  -1.700E+01

 COEFICIENTES DE MARKOV
  -0.5000000    0.5000000
   0.3333333   -0.5000000    0.1666667
  -0.2500000    0.4583333   -0.2500000    0.0416667

 --- RESULTADO FINAL ---
 POTENCIA,    COEFICIENTES
    0         1.0000
    1        -7.0000
    2        15.0000
    3       -10.0000
    4         1.0000
```

✅ **Coincide exactamente con el libro.** $f(1)$ da $1.7\times10^{-16}$ en R contra
$-8.9\times10^{-8}$ en el libro (precisión simple).

| Verificación | Resultado |
|---|---|
| Coeficientes vs. `eigen()` | diferencia máx. $2.1\times10^{-14}$ |
| $c_0 = \det A$ | $1$ ✓ |
| $-c_3 = $ traza | $10$ ✓ |
| Raíces por Frobenius + QR | $0.283119,\ 0.426022,\ 1,\ 8.290859$ ✓ |

**Ejemplo del tema** $x^5 - 0.2x^4 + 7x^3 + x^2 - 3.5x + 2 = 0$ (Frobenius + QR):

```
 Matriz tras la iteracion QR (triangular por bloques):
        [,1]    [,2]    [,3]    [,4]    [,5]
[1,] -2.4159 -6.2305  2.7533  0.1331 -1.1226
[2,]  2.2463  2.6967 -1.3906 -0.0823  0.5659     <- bloque 2x2: par complejo
[3,]  0.0000  0.0000 -0.9086 -0.4991 -0.1163     <- bloque 1x1: real
[4,]  0.0000  0.0000  0.0000  0.4026 -0.8632
[5,]  0.0000  0.0000  0.0000  0.1426  0.4252     <- bloque 2x2: par complejo
```

Raíces: $-0.9086$, $\;0.1404 \pm 2.7314i$, $\;0.4139 \pm 0.3506i$ ✓ (coincide con `polyroot()`
hasta $10^{-14}$). Se ve en vivo la "triangular superior por bloques" del tema.

## Sensibilidad a $\Delta\lambda$

| $\Delta\lambda$ | error máx. en los coeficientes | coef. de $\lambda^4$ |
|---|---:|---:|
| 10 | $1.5\times10^{-12}$ | 1.000000 |
| 2 | $3.8\times10^{-14}$ | 1.000000 |
| **0.5** | $2.1\times10^{-14}$ | 1.000000 |
| 0.1 | $1.0\times10^{-12}$ | 1.000000 |
| 0.01 | $3.3\times10^{-8}$ | 1.000000 |
| 0.001 | $5.6\times10^{-6}$ | 1.000006 |
| 0.0001 | **1.5** | **2.544** |

> **La misma lección del capítulo 5:** con $\Delta\lambda$ muy pequeño, $\Delta^4 f \approx
> 24\,\Delta\lambda^4$ se calcula restando números casi iguales del orden de $1$, y luego se
> divide entre $\Delta\lambda^4$. Con $\Delta\lambda = 10^{-4}$ el coeficiente principal sale
> $2.54$ en lugar de $1$. El libro lo advierte: *"la precisión de los coeficientes puede
> verse afectada por la entrada para el intervalo de $\lambda$"*.

## Programa 7-2

```
 Paso 1:  G = 3.741657   h = 25.224972   u = (0.00000, 6.74166, 2.00000, 1.00000)
 Paso 2:  G = -0.462910  h = 0.374706    u = (0.00000, 0.00000, -0.80946, -0.30690)

 MATRIZ DE HESSENBERG O TRIDIAGONAL
 4.00000E+00 -3.74166E+00  0.00000E+00  0.00000E+00
-3.74166E+00  5.00000E+00  4.62910E-01  0.00000E+00
 0.00000E+00  4.62910E-01  6.66667E-01 -8.90871E-02
 0.00000E+00  0.00000E+00 -8.90871E-02  3.33333E-01

 VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 1 x 1
 4.00000E+00
 VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 2 x 2
 7.25083E-01  8.27492E+00
 VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 3 x 3
 3.81086E-01  9.94724E-01  8.29086E+00
 RESULTADOS FINALES (VALORES PROPIOS DE LA MATRIZ COMPLETA)
 2.83119E-01  4.26022E-01  1.00000E+00  8.29086E+00
```

✅ La tridiagonal y las submatrices coinciden con el libro.

**En la fila final, R es más preciso que el libro:**

| | Libro | R | Exacto (`eigen`) |
|---|---|---|---|
| $\lambda_1$ | 0.283250 | **0.283119** | 0.283119 |
| $\lambda_3$ | 0.999996 | **1.000000** | 1 |

El libro trabaja en precisión simple y corta la iteración con tolerancia $10^{-6}$. Que
$\lambda = 1$ es exacto lo confirma el Programa 7-1: $f(1) = 0$.

| Verificación | Resultado |
|---|---|
| Traza de $A$ y de $T$ | 10 y 10 (se conserva) |
| $\det A$ y $\det T$ | 1 y 1 (se conserva) |
| vs. `eigen()` de R | diferencia máx. $8.3\times10^{-12}$ |
| $\max\lvert PP - I\rvert$ | $1.1\times10^{-16}$ |

**Los dos programas se confirman entre sí:** las raíces del polinomio del 7-1 son
exactamente los valores propios del 7-2.

---

# PARTE VI — GUION DE LA EXPOSICIÓN

## Lámina 1 — Portada
> "Capítulo 7 de Nakamura: Valores Propios. Vamos a ver dos programas que atacan el mismo
> problema por caminos distintos: el 7-1 construye el polinomio característico por
> interpolación, y el 7-2 calcula los valores propios de una matriz simétrica con
> Householder y bisección."

## Lámina 2 — Motivación
> "Los valores propios aparecen cuando un sistema vibra sin fuerzas externas: una cuerda,
> una estructura, un edificio en un sismo. En este sistema de tres masas, los valores
> propios son $\omega^2$, el cuadrado de las frecuencias naturales, y los vectores propios
> son las formas de vibrar: las tres masas juntas, la del centro quieta, o alternadas."

## Lámina 3 — El problema
> "Buscamos $\lambda$ y $\mathbf{x}\ne 0$ con $A\mathbf{x} = \lambda\mathbf{x}$. Eso es un
> sistema homogéneo $(A - \lambda I)\mathbf{x} = 0$, que solo tiene solución no trivial si
> el determinante es cero. Ese determinante es un polinomio de grado $N$ en $\lambda$: el
> polinomio característico. Dos chequeos útiles: la suma de los valores propios es la
> traza y su producto es el determinante."

## Lámina 4 — El método de interpolación
> "Desarrollar $\det(A - \lambda I)$ simbólicamente es inviable para $N$ grande. La idea del
> Programa 7-1 es otra: el polinomio tiene grado $N$, así que si lo evalúo en $N+1$ puntos
> queda determinado exactamente. Cada evaluación es un determinante numérico, que ya sabemos
> hacer con Gauss del capítulo 6. Y con esos $N+1$ valores interpolo con Newton hacia
> adelante, como en el tema 2."

## Lámina 5 — Tabla de diferencias y Markov
> "Con $\Delta\lambda = 0.5$ calculamos 5 determinantes y la tabla de diferencias. La fórmula
> de Newton da el polinomio en binomiales $\binom{s}{k}$. Los coeficientes de Markov son
> simplemente esos binomiales desarrollados en potencias de $s$. Se suman potencia por
> potencia y al final se divide entre $\Delta\lambda^i$ para volver a $\lambda$."

## Lámina 6 — Programa 7-1: salida
> "El resultado es $\lambda^4 - 10\lambda^3 + 15\lambda^2 - 7\lambda + 1$, igual que el libro.
> Se verifica solo: el término independiente es el determinante, 1; el de $\lambda^3$ es
> menos la traza, 10; y $f(1) = 0$, así que 1 es valor propio."

## Lámina 7 — La trampa de $\Delta\lambda$
> "Pero el libro advierte que la precisión depende de $\Delta\lambda$. Lo medimos: con 0.5 el
> error es $10^{-14}$, pero con $10^{-4}$ el coeficiente principal sale 2.5 en lugar de 1.
> Es la misma cancelación del capítulo 5: la cuarta diferencia es como una cuarta derivada,
> y restar números casi iguales destruye las cifras."

## Lámina 8 — Del polinomio a las raíces
> "El 7-1 da el polinomio, no las raíces. Para eso el tema propone la matriz de Frobenius:
> una matriz cuyo polinomio característico es el que queremos. Le aplicamos la iteración QR,
> que la lleva a una triangular por bloques: los bloques de $1\times1$ son raíces reales y
> los de $2\times2$ son pares complejos. Con el polinomio de grado 5 del tema se ve
> clarito."

## Lámina 9 — Householder: la idea
> "El Programa 7-2 va por otro camino, solo para matrices simétricas. Primero la
> tridiagonaliza con transformaciones de Householder, $PAP$ con
> $P = I - \mathbf{u}\mathbf{u}^T/h$. $P$ es simétrica y es su propia inversa, así que es
> una similaridad: no cambia los valores propios. Cada paso anula una columna debajo de la
> subdiagonal; con $N-2$ pasos queda tridiagonal."

## Lámina 10 — Householder a mano
> "En el ejemplo del tema basta un paso. $G = \sqrt{17}$, $h = 21.12$,
> $\mathbf{u} = (0, 5.12, 4)$. El signo de $G$ se elige igual al de $a_{21}$ para que
> $a + G$ sea una suma y no una resta: de nuevo, evitar la cancelación. Y $PP = I$ sale de
> que $\mathbf{u}^T\mathbf{u} = 2h$."

## Lámina 11 — ¿Por qué tridiagonal?
> "Porque el determinante de una tridiagonal sale con una recurrencia de tres términos, en
> $O(N)$ en lugar de $O(N^3)$. Y además $p_i$ es el polinomio característico de la
> submatriz principal $i\times i$, así que de una sola pasada tenemos los de todas."

## Lámina 12 — Entrelazamiento y bisección
> "Lo elegante: los valores propios de la submatriz $L-1$ separan a los de la submatriz $L$.
> Entonces cada intervalo entre dos raíces consecutivas contiene exactamente una raíz nueva.
> El programa empieza por $1\times1$, que es trivial, y va subiendo. En cada intervalo
> aplica bisección mientras es mayor que 1, y después interpolación lineal, como en el
> tema 3."

## Lámina 13 — Programa 7-2: salida
> "La tridiagonal y los valores propios de las submatrices coinciden con el libro. En la
> fila final R es más preciso: el libro da 0.999996 y 0.283250, y los exactos son 1 y
> 0.283119. Es precisión simple con tolerancia $10^{-6}$. Y los dos programas se confirman
> entre sí: las raíces del polinomio del 7-1 son los valores propios del 7-2."

## Lámina 14 — Método de las potencias
> "El tema también presenta el método de las potencias: multiplicar por $A$ una y otra vez.
> El vector se va alineando con el vector propio dominante. Reproducimos la tabla 7.3:
> converge a 8.508. La velocidad la da el cociente $|\lambda_2/\lambda_1|$. Con potencias
> inversas, que es resolver $A\mathbf{v} = \mathbf{u}$ en cada paso, se obtiene el de menor
> valor absoluto."

## Lámina 15 — Comparación
> "Resumiendo los métodos: interpolación sirve para cualquier matriz pero da el polinomio y
> es sensible a $\Delta\lambda$; Householder con bisección es robusto pero solo para
> simétricas; potencias es muy simple pero da solo un valor propio; y QR es el método
> general, el que usan los programas profesionales."

## Lámina 16 — Conclusiones
> "Los valores propios son las raíces del polinomio característico, pero nunca se calculan
> desarrollando el determinante. El 7-1 lo reconstruye por interpolación; el 7-2
> tridiagonaliza con Householder y aprovecha el entrelazamiento para encerrar cada raíz. Y
> aparecen las lecciones de todo el curso: Gauss del capítulo 6, Newton del 2, bisección
> del 3 y la cancelación del 5."

---

# PARTE VII — PREGUNTAS PROBABLES EN LA DEFENSA

**¿Por qué no desarrollar $\det(A - \lambda I)$ directamente?**
Por cofactores cuesta $\sim N!$ operaciones. El 7-1 hace $N+1$ determinantes numéricos de
$O(N^3)$ cada uno.

**¿Por qué la interpolación es exacta?**
Porque $f(\lambda)$ es un polinomio de grado $N$, y un polinomio de grado $N$ queda
determinado de forma única por $N+1$ puntos. El único error es el de redondeo.

**¿Qué son los coeficientes de Markov?**
Los coeficientes de $\binom{s}{k} = \frac{s(s-1)\cdots(s-k+1)}{k!}$ escrito en potencias de
$s$. Sirven para pasar la fórmula de Newton (en binomiales) a un polinomio común.

**¿Cómo elegir $\Delta\lambda$?**
Del orden de los valores propios. Muy pequeño produce cancelación en la tabla de
diferencias; muy grande hace los determinantes enormes. Con esta matriz, entre $0.5$ y $2$
el error es de $10^{-14}$.

**¿Por qué Householder no cambia los valores propios?**
Porque $P = P^{-1}$, así que $PAP = PAP^{-1}$ es una transformación de similaridad:
$\det(PAP^{-1} - \lambda I) = \det(P(A - \lambda I)P^{-1}) = \det(A - \lambda I)$.

**¿Por qué el 7-2 solo sirve para simétricas?**
Porque Householder sobre una simétrica da tridiagonal (la simetría se conserva). Sobre una
no simétrica da Hessenberg, la recurrencia de tres términos no aplica, el entrelazamiento
no está garantizado y los valores propios pueden ser complejos.

**¿Qué garantiza que haya una raíz en cada intervalo?**
El teorema de entrelazamiento para matrices simétricas: los valores propios de la submatriz
principal $L-1$ separan estrictamente a los de la $L$ (si los $b_i \ne 0$).

**¿Por qué bisección y luego interpolación lineal?**
Bisección siempre converge, pero lento (gana un bit por paso). Cuando el intervalo es
pequeño, la interpolación lineal converge más rápido. Se combinan la seguridad de una y la
velocidad de la otra.

**¿Cuándo falla el método de las potencias?**
Si los dos valores propios de mayor módulo son iguales en valor absoluto (por ejemplo
$\lambda$ y $-\lambda$, o un par complejo), o si el vector inicial no tiene componente en
la dirección del vector propio dominante.

**¿Por qué el libro da 0.999996 si el valor exacto es 1?**
Precisión simple (~7 cifras) y tolerancia de $10^{-6}$ en la iteración. En doble precisión
sale 1.000000.

---

# PARTE VIII — RESUMEN COMPARATIVO

| Método | Matriz | Da | Costo | Punto débil |
|---|---|---|---|---|
| **Interpolación (7-1)** | cualquiera | polinomio característico | $(N+1)$ dets $\to O(N^4)$ | sensible a $\Delta\lambda$; falta hallar raíces |
| **Householder + bisección (7-2)** | simétrica | todos los $\lambda$ (reales) | $O(N^3)$ + $O(N)$ por evaluación | solo simétricas |
| Potencias | cualquiera | el de mayor $\lvert\lambda\rvert$ | $O(N^2)$ por iteración | uno solo; lento si $\lvert\lambda_2/\lambda_1\rvert \approx 1$ |
| Potencias inversas | cualquiera | el de menor $\lvert\lambda\rvert$ | resolver un sistema por iteración | uno solo |
| QR | cualquiera | todos (incluso complejos) | $O(N^3)$ | el más complejo de programar |

## Decisión rápida

```
¿La matriz es simétrica?
├── Sí → HOUSEHOLDER + BISECCIÓN (Programa 7-2)
└── No → QR (o interpolación + raíces si N es chico)

¿Necesito solo el mayor (o el menor)?
└── Sí → POTENCIAS (o POTENCIAS INVERSAS)

¿Necesito el polinomio característico en sí?
└── Sí → INTERPOLACIÓN (Programa 7-1), con Δλ del orden de los λ

¿Quiero las raíces de un polinomio?
└── Matriz de FROBENIUS + QR
```

> **Frase para cerrar:** *"Los valores propios son raíces de un polinomio, pero el polinomio
> nunca se escribe a mano: o se reconstruye por interpolación, o se evita del todo
> transformando la matriz."*
