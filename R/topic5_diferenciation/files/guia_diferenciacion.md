# Guía de Diferenciación Numérica — Capítulo 5 (Nakamura)

> Universidad de Los Andes · CeSiMo · Análisis Numérico
> Programa 5-1 implementado en R · Defensa: martes 29 de septiembre

---

# PARTE I — FORMULARIO

## 1. Primeras derivadas — fórmulas básicas

| Nombre | Fórmula | Error |
|---|---|:---:|
| **Adelante** (2 pts) | $f'_i \approx \dfrac{f_{i+1} - f_i}{h}$ | $O(h)$ |
| **Atrás** (2 pts) | $f'_i \approx \dfrac{f_i - f_{i-1}}{h}$ | $O(h)$ |
| **Central** (2 pts) | $f'_i \approx \dfrac{f_{i+1} - f_{i-1}}{2h}$ | $O(h^2)$ |
| **Adelante** (3 pts) | $f'_i \approx \dfrac{-3f_i + 4f_{i+1} - f_{i+2}}{2h}$ | $O(h^2)$ |
| **Atrás** (3 pts) | $f'_i \approx \dfrac{3f_i - 4f_{i-1} + f_{i-2}}{2h}$ | $O(h^2)$ |
| **Central** (5 pts) | $f'_i \approx \dfrac{f_{i-2} - 8f_{i-1} + 8f_{i+1} - f_{i+2}}{12h}$ | $O(h^4)$ |

**Términos de error explícitos:**

$$
\begin{aligned}
\text{Adelante 2 pts:}&\quad E = -\tfrac{h}{2}f''_i \\
\text{Atrás 2 pts:}&\quad E = +\tfrac{h}{2}f''_i \\
\text{Central 2 pts:}&\quad E = -\tfrac{h^2}{6}f'''_i \\
\text{Adelante 3 pts:}&\quad E = -\tfrac{h^2}{3}f'''_i \\
\text{Atrás 3 pts:}&\quad E = +\tfrac{h^2}{3}f'''_i
\end{aligned}
$$

---

## 2. Segundas derivadas

| Nombre | Fórmula | Error |
|---|---|:---:|
| **Central** (3 pts) | $f''_i \approx \dfrac{f_{i+1} - 2f_i + f_{i-1}}{h^2}$ | $O(h^2)$ |
| **Atrás** (3 pts) | $f''_i \approx \dfrac{f_{i-2} - 2f_{i-1} + f_i}{h^2}$ | $O(h)$ |
| **Central** (5 pts) | $f''_i \approx \dfrac{-f_{i-2} + 16f_{i-1} - 30f_i + 16f_{i+1} - f_{i+2}}{12h^2}$ | $O(h^4)$ |

Errores: central 3 pts → $-\tfrac{h^2}{12}f^{(4)}_i$ &nbsp;·&nbsp; atrás 3 pts → $h\,f'''_i$

---

## 3. Regla mínima de puntos

$$\boxed{\;\text{Para una derivada de orden } p \text{ se necesitan al menos } p+1 \text{ puntos}\;}$$

Con $L$ puntos se pueden imponer $L$ condiciones, así que el error queda
proporcional al término $(L+1)$-ésimo del desarrollo de Taylor.

---

## 4. Fórmula general (ecuación 5.3.1)

$$\boxed{\;f^{(p)}_0 = \frac{a_\alpha f_\alpha + a_\beta f_\beta + \cdots + a_\lambda f_\lambda}{h^p} + E\;}$$

con el término de error (ecuación 5.3.2):

$$\boxed{\;E = c_1\,h^{L-p}f^{(L)} + c_2\,h^{L-p+1}f^{(L+1)}\;}$$

donde $x_i = ih$ para $i = \alpha, \beta, \dots, \lambda$, y la derivada se
evalúa en $x = 0$.

---

## 5. Sistema lineal que determina los coeficientes

Sustituyendo los desarrollos de Taylor y agrupando por derivada, se obtiene
un sistema de **Vandermonde**:

$$
\begin{bmatrix}
1 & 1 & \cdots & 1 \\
\alpha & \beta & \cdots & \lambda \\
\alpha^2 & \beta^2 & \cdots & \lambda^2 \\
\vdots & \vdots & & \vdots \\
\alpha^{L-1} & \beta^{L-1} & \cdots & \lambda^{L-1}
\end{bmatrix}
\begin{bmatrix} a_\alpha \\ a_\beta \\ \vdots \\ a_\lambda \end{bmatrix}
=
\begin{bmatrix} 0 \\ \vdots \\ p! \\ \vdots \\ 0 \end{bmatrix}
\;\leftarrow\; \text{fila } p+1
$$

**En palabras:** la fila $k$ contiene las potencias $(k-1)$ de los índices,
y el lado derecho es todo ceros salvo $p!$ en la fila $p+1$.

**Coeficientes del error:**

$$c_1 = -\frac{1}{L!}\sum_j a_j\,j^{\,L}, \qquad
c_2 = -\frac{1}{(L+1)!}\sum_j a_j\,j^{\,L+1}$$

> Si $c_1 = 0$, el error es del orden siguiente y manda $c_2$.

---

## 6. Operadores de diferencias

| Operador | Símbolo | Definición |
|---|:---:|---|
| Hacia adelante | $\Delta$ | $\Delta f_i = f_{i+1} - f_i$ |
| Hacia atrás | $\nabla$ | $\nabla f_i = f_i - f_{i-1}$ |
| Central | $\delta$ | $\delta f_i = f_{i+1/2} - f_{i-1/2}$ |

**Orden 2:**

$$\Delta^2 f_i = f_{i+2} - 2f_{i+1} + f_i$$
$$\nabla^2 f_i = f_i - 2f_{i-1} + f_{i-2}$$
$$\delta^2 f_i = f_{i+1} - 2f_i + f_{i-1}$$

**Identidad clave:**

$$\boxed{\;\delta^2 = \Delta\nabla = \nabla\Delta\;}$$

Aproximación de operadores diferenciales:
$\dfrac{d}{dx} \approx \dfrac{\Delta}{\Delta x} \approx \dfrac{\nabla}{\nabla x} \approx \dfrac{\delta}{\delta x}$

---

## 7. Derivadas parciales

**Primeras (central):**

$$f_x \approx \frac{f(x_0+\Delta x,\,y_0) - f(x_0-\Delta x,\,y_0)}{2\Delta x}$$

**Segundas (central):**

$$f_{xx} \approx \frac{f(x_0+\Delta x,y_0) - 2f(x_0,y_0) + f(x_0-\Delta x,y_0)}{(\Delta x)^2}$$

$$f_{yy} \approx \frac{f(x_0,y_0+\Delta y) - 2f(x_0,y_0) + f(x_0,y_0-\Delta y)}{(\Delta y)^2}$$

**Cruzada:**

$$f_{xy} \approx \frac{f(x_0{+}\Delta x,y_0{+}\Delta y) - f(x_0{+}\Delta x,y_0{-}\Delta y) - f(x_0{-}\Delta x,y_0{+}\Delta y) + f(x_0{-}\Delta x,y_0{-}\Delta y)}{4\Delta x\,\Delta y}$$

> Regla mnemotécnica de $f_{xy}$: los signos van $+\;-\;-\;+$ siguiendo las
> cuatro esquinas del cuadrado, como un determinante $2\times2$.

---

# PARTE II — TEORÍA

## De dónde salen las fórmulas: el desarrollo de Taylor

Todas las aproximaciones por diferencias nacen del mismo lugar. El desarrollo
de Taylor alrededor de $x_i$:

$$f_{i+1} = f_i + hf'_i + \frac{h^2}{2}f''_i + \frac{h^3}{6}f'''_i + \frac{h^4}{24}f^{(4)}_i + \cdots$$

### Diferencia hacia adelante

Despejando $f'_i$ directamente:

$$f'_i = \frac{f_{i+1} - f_i}{h} - \underbrace{\frac{h}{2}f''_i - \frac{h^2}{6}f'''_i - \cdots}_{\text{error}}$$

$$\Rightarrow\quad f'_i = \frac{f_{i+1} - f_i}{h} + O(h)$$

El primer término descartado es $-\frac{h}{2}f''_i$ → error de **primer orden**.

### Diferencia central: por qué gana un orden

Escribimos también hacia atrás:

$$f_{i-1} = f_i - hf'_i + \frac{h^2}{2}f''_i - \frac{h^3}{6}f'''_i + \cdots$$

Al **restar** $f_{i+1} - f_{i-1}$, los términos pares se cancelan:

$$f_{i+1} - f_{i-1} = 2hf'_i + \frac{h^3}{3}f'''_i + \cdots$$

$$\Rightarrow\quad f'_i = \frac{f_{i+1} - f_{i-1}}{2h} - \frac{h^2}{6}f'''_i$$

> **La clave:** el término en $h^2$ (que contenía $f''$) **desaparece por
> simetría**. Por eso la diferencia central es $O(h^2)$ con los mismos dos
> puntos que las laterales usan para lograr solo $O(h)$.

### Segunda derivada central

Ahora **sumamos** en vez de restar — se cancelan los términos impares:

$$f_{i+1} + f_{i-1} = 2f_i + h^2f''_i + \frac{h^4}{12}f^{(4)}_i + \cdots$$

$$\Rightarrow\quad f''_i = \frac{f_{i+1} - 2f_i + f_{i-1}}{h^2} - \frac{h^2}{12}f^{(4)}_i$$

### Adelante con 3 puntos: eliminar términos a mano

Con dos desarrollos:

$$f_{i+1} = f_i + hf'_i + \tfrac{h^2}{2}f''_i + \tfrac{h^3}{6}f'''_i + \cdots$$
$$f_{i+2} = f_i + 2hf'_i + \tfrac{4h^2}{2}f''_i + \tfrac{8h^3}{6}f'''_i + \cdots$$

Combinando $4f_{i+1} - f_{i+2}$ el término de $f''$ se anula:

$$4f_{i+1} - f_{i+2} = 3f_i + 2hf'_i - \tfrac{2}{3}h^3f'''_i + \cdots$$

$$\Rightarrow\quad f'_i = \frac{-f_{i+2} + 4f_{i+1} - 3f_i}{2h} - \frac{h^2}{3}f'''_i$$

Se logra $O(h^2)$ **sin usar puntos a la izquierda** — útil en un borde donde
no hay datos hacia atrás.

---

## El algoritmo genérico (sección 5.3)

Los ejemplos anteriores se resolvieron "a mano", combinando desarrollos hasta
que se cancelara lo que estorbaba. El **Programa 5-1** automatiza eso para
cualquier configuración.

### Planteo

Se quiere:

$$f^{(p)}_0 = \frac{a_\alpha f_\alpha + \cdots + a_\lambda f_\lambda}{h^p} + E$$

Se sustituye cada $f_j$ por su desarrollo de Taylor alrededor de $x=0$ y se
agrupa por derivada.

### Ejemplo desarrollado: $p=1$, $L=3$, puntos $0,1,2$

$$f'_0 = \frac{1}{h}\Bigl[a_0 f_0 + a_1\bigl(f_0 + hf'_0 + \tfrac{h^2}{2}f''_0 + \cdots\bigr) + a_2\bigl(f_0 + 2hf'_0 + \tfrac{4h^2}{2}f''_0 + \cdots\bigr)\Bigr] + E$$

Reagrupando por derivada:

$$f'_0 = \frac{f_0[a_0+a_1+a_2]}{h} + f'_0[a_1 + 2a_2] + \frac{h\,f''_0[a_1+4a_2]}{2} + \frac{h^2 f'''_0[a_1+8a_2]}{6} + \cdots$$

Para que la igualdad se cumpla necesitamos que el coeficiente de $f_0$ sea $0$,
el de $f'_0$ sea $1$, y el de $f''_0$ sea $0$:

$$
\begin{cases}
a_0 + a_1 + a_2 = 0 & \text{(coef. de } f_0)\\
\phantom{a_0 +\;} a_1 + 2a_2 = 1 & \text{(coef. de } f'_0)\\
\phantom{a_0 +\;} a_1 + 4a_2 = 0 & \text{(coef. de } f''_0)
\end{cases}
$$

**Solución:** $a_0 = -\tfrac32$, $a_1 = 2$, $a_2 = -\tfrac12$

$$\Rightarrow\quad f'_0 = \frac{-3f_0 + 4f_1 - f_2}{2h} + E$$

### El error sale de las filas que sobran

Los términos que **no** anulamos forman el error:

$$c_1 = -\frac{1}{6}(a_1 + 8a_2) = -\frac{1}{6}(2 - 4) = \frac{1}{3}$$
$$c_2 = -\frac{1}{24}(a_1 + 16a_2) = -\frac{1}{24}(2 - 8) = \frac{1}{4}$$

$$E = \tfrac13 h^2 f'''_0 + \tfrac14 h^3 f^{(4)}_0 + \cdots$$

### Generalización

Con $L$ puntos se fijan los primeros $L$ términos del desarrollo de Taylor.
El error queda proporcional al término $(L+1)$-ésimo — o al siguiente si el
coeficiente se anula.

**Ventaja enorme:** el algoritmo funciona aunque los índices $\alpha,\beta,\dots$
**no sean enteros**. Eso significa que sirve para retículas con espaciamiento
**no uniforme**.

---

# PARTE III — EL PROGRAMA 5-1 PASO A PASO

## Qué hace (y qué no hace)

> **No calcula una derivada.** Recibe una configuración de puntos y un orden $p$,
> y **devuelve la fórmula** de diferencias correspondiente junto con su término
> de error. Es un *generador de fórmulas*.

**Entradas:**
1. $L$ = número de puntos de la retícula (máx. 10)
2. Los índices $\alpha, \beta, \dots, \lambda$
3. $p$ = orden de la derivada

**Salidas:**
1. Los coeficientes $a_\alpha, \dots, a_\lambda$ en forma racional
2. Los coeficientes $c_1$ y $c_2$ del término de error

---

## Estructura del código R

```
gauss_solve()   → eliminación gaussiana con pivoteo (equiv. SUBROUTINE GAUSS)
to_frac()       → convierte decimales a fracciones (fracción continua)
diff_scheme()   → construye el sistema, lo resuelve, calcula el error
print_scheme()  → salida en el formato del libro
print_fraction()→ forma racional exacta
verify()        → comprueba la fórmula contra una derivada conocida
```

---

## `diff_scheme()` — el corazón

### Paso 1: construir la matriz de Vandermonde

```r
B <- outer(0:(L + 1), el, function(k, e) e^k)
B[1, ] <- 1          # fila k=1 -> el^0 = 1 (evita el caso 0^0)
```

Esto genera $L+2$ filas: las primeras $L$ son el sistema, y las **dos últimas
se guardan para calcular el error**.

Para `el = c(0,1,2)`:

| fila $k$ | potencia | $j=0$ | $j=1$ | $j=2$ |
|:---:|:---:|:---:|:---:|:---:|
| 1 | $j^0$ | 1 | 1 | 1 |
| 2 | $j^1$ | 0 | 1 | 2 |
| 3 | $j^2$ | 0 | 1 | 4 |
| *4* | *$j^3$* | *0* | *1* | *8* | ← error
| *5* | *$j^4$* | *0* | *1* | *16* | ← error

> **En el FORTRAN esto es:** `A(K,L)=EL(L)**(K-1)` guardado también en `B(K,L)`
> porque `A` se destruye al resolver el sistema y `B` se necesita después.

### Paso 2: lado derecho

```r
rhs <- numeric(L)
rhs[p + 1] <- factorial(p)
```

Todo ceros salvo $p!$ en la posición $p+1$. Con $p=1$: `rhs = c(0, 1, 0)`.

> **En el FORTRAN:** `IF (K-1.EQ.KDR) A(K,KM+1)=Z` donde `Z` se calculó como
> $p!$ en el bucle `DO I=1,KDR`.

### Paso 3: resolver

```r
a <- gauss_solve(A_sys, rhs)
```

### Paso 4: coeficientes del error

```r
C1 <- sum(B[L + 1, ] * a)      # fila L+1
C2 <- sum(B[L + 2, ] * a)      # fila L+2
c1 <- -C1 / factorial(L)
c2 <- -C2 / factorial(L + 1)
```

> **En el FORTRAN:** el bucle `C(K)=C(K)+B(K,L)*A(L,KM+1)` — usa `B` (la copia)
> multiplicada por la solución que quedó en la última columna de `A`.

---

## `gauss_solve()` — eliminación gaussiana

Se implementó desde cero para ser fiel al listado (el libro la explica recién
en el capítulo 6).

```r
# --- Eliminación hacia adelante con pivoteo parcial ---
for (i in seq_len(n - 1)) {
  ipv <- i
  for (j in (i + 1):n)
    if (abs(M[ipv, i]) < abs(M[j, i])) ipv <- j    # buscar mayor pivote
  if (ipv != i) M[c(i, ipv), ] <- M[c(ipv, i), ]   # intercambiar filas

  for (jr in (i + 1):n) {
    r <- M[jr, i] / M[i, i]
    M[jr, i:(n+1)] <- M[jr, i:(n+1)] - r * M[i, i:(n+1)]
  }
}

# --- Sustitución hacia atrás ---
x[n] <- M[n, n+1] / M[n, n]
for (nv in (n-1):1)
  x[nv] <- (M[nv, n+1] - sum(M[nv, (nv+1):n] * x[(nv+1):n])) / M[nv, nv]
```

**¿Por qué pivoteo parcial?** Porque la matriz de Vandermonde tiene ceros en la
primera columna cuando algún índice es $0$ (ej. `el = c(0,1,2)` → columna
$[1,0,0]^T$). Sin pivoteo se dividiría entre cero.

---

## La normalización racional

Los coeficientes salen como decimales ($-1.5,\, 2,\, -0.5$). Para mostrarlos
como enteros, el libro divide por el menor valor absoluto no nulo:

```r
F  <- min(abs(a[abs(a) > 1e-4]))    # F = 0.5
CF <- a / F                          # [-3, 4, -1]
# denominador mostrado = 1/F = 2
```

En R se agregó además una versión **exacta** con fracciones continuas, que
funciona incluso cuando la normalización del libro no da enteros:

```r
fr  <- lapply(a, to_frac)                       # cada a_j como num/den
den <- Reduce(lcm2, sapply(fr, `[`, 2))         # mínimo común denominador
num <- round(a * den)                           # numeradores enteros
```

---

# PARTE IV — SALIDAS VERIFICADAS

## Ejemplo 1 del libro: 3 puntos $(0,1,2)$, $p=1$

```
 ESQUEMA DE DIFERENCIAS
 +[  -3.00000/( 2.00000 H**1)] F( 0.000H)
 +[   4.00000/( 2.00000 H**1)] F( 1.000H)
 +[  -1.00000/( 2.00000 H**1)] F( 2.000H)

 TERMINO DEL ERROR
     (   2.00000/   6.00000)H**2  F^(3)
     (   6.00000/  24.00000)H**3  F^(4)
```

✅ **Coincide exactamente con el libro.**

$$f'_0 = \frac{-3f_0 + 4f_1 - f_2}{2h} + \tfrac13 h^2 f'''$$

**Verificación numérica** con $f(x) = e^x$ (donde $f'(0) = 1$):

| $h$ | aproximado | error | ratio |
|---|---|---|:---:|
| 0.1000 | 0.99640457 | 3.60e−03 | — |
| 0.0500 | 0.99913467 | 8.65e−04 | 4.2 |
| 0.0250 | 0.99978771 | 2.12e−04 | 4.1 |

El error se divide por **4** al partir $h$ a la mitad → confirma $O(h^2)$.

---

## Ejemplo 2 del libro: 5 puntos $(-2..2)$, $p=2$

```
 ESQUEMA DE DIFERENCIAS
 +[  -1.00000/(12.00000 H**2)] F(-2.000H)
 +[  16.00000/(12.00000 H**2)] F(-1.000H)
 +[ -30.00000/(12.00000 H**2)] F( 0.000H)
 +[  16.00000/(12.00000 H**2)] F( 1.000H)
 +[  -1.00000/(12.00000 H**2)] F( 2.000H)

 TERMINO DEL ERROR
     (   0.00000/ 120.00000)H**3  F^(5)
     (   8.00000/ 720.00000)H**4  F^(6)
```

✅ **Coincide exactamente con el libro.**

Nótese que $c_1 = 0$: el término en $h^3$ **se anula por simetría**, así que
el error real es $O(h^4)$, no $O(h^3)$.

| $h$ | aproximado | error | ratio |
|---|---|---|:---:|
| 0.1000 | 0.99999889 | 1.11e−06 | — |
| 0.0500 | 0.99999993 | 6.95e−08 | 16.0 |
| 0.0250 | 1.00000000 | 4.34e−09 | 16.0 |

El error se divide por **16** → confirma $O(h^4)$.

---

## Fórmulas clásicas regeneradas

El mismo algoritmo produce todas las fórmulas de la tabla 5.2:

| Configuración | Fórmula generada | Orden |
|---|---|:---:|
| $(0,1)$, $p{=}1$ | $[-f_0 + f_1]/h$ | $O(h)$ |
| $(-1,0)$, $p{=}1$ | $[-f_{-1} + f_0]/h$ | $O(h)$ |
| $(-1,1)$, $p{=}1$ | $[-f_{-1} + f_1]/2h$ | $O(h^2)$ |
| $(0,1,2)$, $p{=}1$ | $[-3f_0 + 4f_1 - f_2]/2h$ | $O(h^2)$ |
| $(-2,-1,0)$, $p{=}1$ | $[f_{-2} - 4f_{-1} + 3f_0]/2h$ | $O(h^2)$ |
| $(-1,0,1)$, $p{=}2$ | $[f_{-1} - 2f_0 + f_1]/h^2$ | $O(h^2)$ |
| $(-2,-1,0)$, $p{=}2$ | $[f_{-2} - 2f_{-1} + f_0]/h^2$ | $O(h)$ |
| $(-2,-1,0,1,2)$, $p{=}1$ | $[f_{-2} - 8f_{-1} + 8f_1 - f_2]/12h$ | $O(h^4)$ |

> **Observación para la exposición:** las configuraciones **simétricas** ganan
> un orden respecto a las asimétricas con el mismo número de puntos. Comparar
> las filas 6 y 7: ambas usan 3 puntos para $f''$, pero la central da $O(h^2)$
> y la de atrás solo $O(h)$.

---

## Retícula NO uniforme

El caso del enunciado: $f''(0)$ con puntos en $-2$, $0.5$ y $1.5$.

```
 Forma racional:  f^(2)(0) = [ +8*f(-2h) -28*f(0.5h) +20*f(1.5h) ] / (35 h^2)
 Error dominante: E = (-13/48) h^2 f^(4)   ->  O(h^2)
```

Ningún libro trae esta fórmula en tablas — el algoritmo la genera igual.
Es la razón de ser del Programa 5-1.

---

# PARTE V — EJEMPLO 5.2 RESUELTO (derivadas parciales)

Tabla de $f(x,y)$ con $\Delta x = \Delta y = h = 0.5$:

| $y \backslash x$ | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 |
|---|---|---|---|---|---|
| **1.0** | 1.63 | 2.05 | 2.50 | 2.98 | 3.49 |
| **1.5** | 1.98 | 2.51 | 3.08 | 3.69 | 4.33 |
| **2.0** | 2.28 | 2.91 | **3.61** | 4.37 | 5.17 |
| **2.5** | 2.64 | 3.25 | 4.08 | 5.00 | 5.98 |
| **3.0** | 2.65 | 3.50 | 4.48 | 5.57 | 6.76 |

## a) Diferencias centrales en $(2,2)$

$$f_x(2,2) = \frac{f(2.5,\,2) - f(1.5,\,2)}{2h} = \frac{4.37 - 2.91}{1.0} = \mathbf{1.46}$$

$$f_y(2,2) = \frac{f(2,\,2.5) - f(2,\,1.5)}{2h} = \frac{4.08 - 3.08}{1.0} = \mathbf{1.00}$$

$$f_{yy}(2,2) = \frac{f(2,2.5) - 2f(2,2) + f(2,1.5)}{h^2} = \frac{4.08 - 2(3.61) + 3.08}{0.25} = \mathbf{-0.24}$$

$$f_{xy}(2,2) = \frac{5.00 - 3.69 - 3.25 + 2.51}{(2\cdot 0.5)^2} = \frac{0.57}{1.0} = \mathbf{0.57}$$

> Para $f_{xy}$ los cuatro valores son las **esquinas** del cuadrado alrededor
> de $(2,2)$, con signos $+\,-\,-\,+$.

## b) Diferencias hacia adelante de 3 puntos

$$f_x(2,2) = \frac{-f(3,2) + 4f(2.5,2) - 3f(2,2)}{2h} = \frac{-5.17 + 4(4.37) - 3(3.61)}{1.0} = \mathbf{1.48}$$

$$f_y(2,2) = \frac{-f(2,3) + 4f(2,2.5) - 3f(2,2)}{2h} = \frac{-4.48 + 4(4.08) - 3(3.61)}{1.0} = \mathbf{1.01}$$

Los resultados de (a) y (b) son muy parecidos ($1.46$ vs $1.48$) porque **ambas
fórmulas son $O(h^2)$**.

---

# PARTE VI — GUION DE LA EXPOSICIÓN

## Lámina 1 — Portada
> "Capítulo 5 de Nakamura: Diferenciación Numérica. Vamos a ver cómo calcular
> derivadas cuando no tenemos la función analítica, y el Programa 5-1, que es
> distinto a todos los anteriores porque **genera fórmulas** en vez de números."

## Lámina 2 — Motivación
> "El problema: tenemos $f$ solo en puntos discretos — una tabla de datos de un
> sensor, o el resultado de una simulación. No podemos derivar analíticamente.
> La idea es aproximar la derivada con combinaciones de los valores conocidos.
> Y a diferencia de la integración, acá hay un peligro nuevo: la diferenciación
> **amplifica el ruido**."

## Lámina 3 — Las tres diferencias básicas
> "Adelante, atrás y central. Las dos primeras son $O(h)$; la central es
> $O(h^2)$ **con los mismos dos puntos**. Eso no es casualidad y lo vemos en la
> lámina siguiente."

## Lámina 4 — Taylor: de dónde salen
> "Todo sale del desarrollo de Taylor. Si despejo directo me queda la diferencia
> hacia adelante, y el primer término que tiro es $-\frac{h}{2}f''$ — de ahí el
> $O(h)$. Ahora, si escribo también hacia atrás y **resto** las dos series, los
> términos pares se cancelan y me queda $O(h^2)$. La simetría es lo que da la
> precisión extra. Y si en vez de restar **sumo**, se cancelan los impares y me
> queda la fórmula de la segunda derivada."

## Lámina 5 — Tabla de fórmulas
> "Estas son las fórmulas de la tabla 5.2 del libro. Lo importante es notar el
> patrón: las configuraciones **simétricas** ganan un orden. Comparen las dos
> fórmulas de $f''$ con 3 puntos: la central da $O(h^2)$, la de atrás solo
> $O(h)$."

## Lámina 6 — El algoritmo genérico
> "Hasta acá derivamos las fórmulas a mano, combinando series de Taylor hasta
> que se cancelara lo que estorbaba. La sección 5.3 automatiza eso. Se plantea
> la fórmula con coeficientes indeterminados, se sustituyen los desarrollos de
> Taylor, se agrupa por derivada, y se pide que el coeficiente de $f^{(p)}$ sea
> 1 y todos los demás 0. Eso da un **sistema lineal**."

## Lámina 7 — El sistema de Vandermonde
> "El sistema tiene una estructura muy limpia: la fila $k$ son las potencias
> $k-1$ de los índices, y el lado derecho es todo ceros salvo $p!$ en la fila
> $p+1$. Es una matriz de **Vandermonde**. Con $L$ puntos tengo $L$ ecuaciones
> y $L$ incógnitas."

## Lámina 8 — El Programa 5-1
> "Acá está la implementación. Dos detalles: primero, se construyen $L+2$ filas
> aunque el sistema use solo $L$ — las dos extra son para calcular el término
> del error. Segundo, se guarda una copia de la matriz en `B` porque la
> eliminación gaussiana destruye `A`. Eso está literal en el FORTRAN del libro."

## Lámina 9 — Eliminación gaussiana
> "El sistema se resuelve con eliminación gaussiana con **pivoteo parcial**.
> El pivoteo acá no es opcional: cuando uno de los índices es 0, la primera
> columna queda con ceros y sin pivoteo se dividiría entre cero. Esto lo
> veremos en detalle en el capítulo 6."

## Lámina 10 — Salidas verificadas
> "Los dos ejemplos del libro reproducidos exactamente. Y además una
> verificación numérica: aplicando la fórmula a $e^x$, el error se divide por 4
> al partir $h$ a la mitad en el primer caso, y por 16 en el segundo —
> confirmando $O(h^2)$ y $O(h^4)$ respectivamente."

## Lámina 11 — Retícula no uniforme
> "Acá está la razón de ser del programa. Con puntos en $-2$, $0.5$ y $1.5$ —
> espaciamiento irregular — genera una fórmula que no está en ninguna tabla:
> $[8f_{-2} - 28f_{0.5} + 20f_{1.5}]/35h^2$. Ese es el valor real de tener un
> generador en vez de una tabla."

## Lámina 12 — Operadores de diferencias
> "Una notación compacta: $\Delta$ hacia adelante, $\nabla$ hacia atrás,
> $\delta$ central. Lo elegante es que los operadores de orden superior son
> **potencias** de estos, y se cumple la identidad $\delta^2 = \Delta\nabla =
> \nabla\Delta$ — las tres formas de escribir la segunda diferencia dan lo
> mismo."

## Lámina 13 — Derivadas parciales
> "La extensión a varias variables es directa: se aplica la misma fórmula en
> cada dirección. Lo único nuevo es la derivada cruzada $f_{xy}$, que usa las
> cuatro esquinas del cuadrado con signos $+\,-\,-\,+$."

## Lámina 14 — El peligro: cancelación
> "Y acá el punto más importante y contraintuitivo. En integración, achicar $h$
> siempre mejora. En diferenciación **no**: al restar dos números casi iguales
> se pierden cifras significativas. Hay un $h$ óptimo, y por debajo de él el
> error **crece**. Es exactamente lo contrario a la intuición."

## Lámina 15 — Conclusiones
> "Resumiendo: las fórmulas salen de Taylor; la simetría regala un orden de
> precisión; el Programa 5-1 automatiza la deducción resolviendo un sistema de
> Vandermonde; funciona incluso con retículas no uniformes; y hay que cuidar
> el $h$ porque muy chico es tan malo como muy grande."

---

# PARTE VII — ADVERTENCIA: EL PELIGRO DE $h$ PEQUEÑO

Este es el punto más importante del capítulo y el que más se pregunta.

## En integración, $h$ pequeño siempre es mejor. En diferenciación, NO.

$$f'(x) \approx \frac{f(x+h) - f(x)}{h}$$

Cuando $h \to 0$, los valores $f(x+h)$ y $f(x)$ se vuelven **casi idénticos**.
Al restarlos se pierden cifras significativas (*cancelación catastrófica*), y
después se divide por un $h$ minúsculo, que **amplifica** el ruido.

## Los dos errores compiten

| Fuente | Comportamiento | Efecto al achicar $h$ |
|---|---|---|
| **Truncamiento** | $\propto h$ (Taylor) | ↓ mejora |
| **Redondeo** | $\propto \varepsilon/h$ | ↑ empeora |

$$E_{\text{total}} \approx \underbrace{\frac{h}{2}|f''|}_{\text{truncamiento}} + \underbrace{\frac{2\varepsilon}{h}}_{\text{redondeo}}$$

Minimizando respecto a $h$:

$$\boxed{\;h_{\text{óptimo}} \approx \sqrt{\varepsilon} \approx 10^{-8}\;}$$

para doble precisión ($\varepsilon \approx 2.2\times10^{-16}$).

Para la diferencia **central**, el óptimo es $h \approx \varepsilon^{1/3} \approx 6\times10^{-6}$.

## Contraste directo

```
INTEGRACIÓN                    DIFERENCIACIÓN
h ↓  →  error ↓  siempre       h ↓  →  error ↓ ... y después ↑
(sumar valores es estable)     (restar valores casi iguales es inestable)
```

> **Frase para la defensa:** *"La integración promedia y por eso suaviza el
> ruido; la diferenciación resta y por eso lo amplifica."*

---

# PARTE VIII — RESUMEN COMPARATIVO

| Fórmula | Puntos | Orden | ¿Cuándo usarla? |
|---|:---:|:---:|---|
| Adelante 2 pts | 2 | $O(h)$ | Borde izquierdo, sin datos previos |
| Atrás 2 pts | 2 | $O(h)$ | Borde derecho, sin datos posteriores |
| **Central 2 pts** | 2 | $O(h^2)$ | **Por defecto** en puntos interiores |
| Adelante 3 pts | 3 | $O(h^2)$ | Borde, pero con precisión $O(h^2)$ |
| Central 5 pts | 5 | $O(h^4)$ | Máxima precisión, función suave |
| Central $f''$ 3 pts | 3 | $O(h^2)$ | Segunda derivada estándar |

## Decisión rápida

```
¿Dónde está el punto?
├── Interior           → CENTRAL (mejor orden con menos puntos)
├── Borde izquierdo    → ADELANTE 3 pts (O(h²) sin mirar a la izquierda)
└── Borde derecho      → ATRÁS 3 pts

¿La retícula es uniforme?
├── Sí   → usar tabla 5.2
└── No   → PROGRAMA 5-1 (genera la fórmula que haga falta)

¿Los datos tienen ruido?
└── Sí   → NO achicar h. Suavizar primero (splines) y derivar el spline.
```
