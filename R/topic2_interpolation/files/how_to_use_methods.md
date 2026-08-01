# **Guía de uso — Interpolación Polinomial**
## <font color="#00b0f0">1. Lagrange</font>
**Cuándo usarlo**: cualquier conjunto de puntos, sin requisitos de espaciado. Ideal cuando tienes pocos puntos (3-5) y necesitas evaluar en un punto específico.

**Pasos**:
1. Para cada $k$, calcular $L_k(x)$: multiplicar los factores $\frac{x - x_j}{x_k - x_j}$ para todo $j \neq k$.
2. Multiplicar $L_k(x) \cdot f(x_k)$.
3. Sumar todos los términos.

**Ejemplo** con 3 puntos: $(1, 0.671),\ (2, 0.620),\ (3, 0.567)$, evaluar en $x = 1.5$
$$L_0(1.5) = \frac{(1.5-2)(1.5-3)}{(1-2)(1-3)} = \frac{(-0.5)(-1.5)}{(-1)(-2)} = \frac{0.75}{2} = 0.375$$
$$L_1(1.5) = \frac{(1.5-1)(1.5-3)}{(2-1)(2-3)} = \frac{(0.5)(-1.5)}{(1)(-1)} = 0.75$$
$$L_2(1.5) = \frac{(1.5-1)(1.5-2)}{(3-1)(3-2)} = \frac{(0.5)(-0.5)}{(2)(1)} = -0.125$$
$$P(1.5) = 0.375(0.671) + 0.75(0.620) + (-0.125)(0.567) = 0.6463$$

**Verificación rápida**: $L_0 + L_1 + L_2 = 1$ siempre. Si no suma 1, hay un error.

---
## <font color="#00b0f0">2. Newton hacia adelante</font>
**Cuándo usarlo**: datos **equiespaciados** ($h$ constante), y el punto $x$ que quieres está cerca del **inicio** de la tabla ($x_0$).

**Pasos**:
1. Construir la tabla de diferencias finitas $\Delta$:

| $i$ | $x_i$ | $f_i$ | $\Delta f$ | $\Delta^2 f$              | $\Delta^3 f$ |
| --- | ----- | ----- | ---------- | ------------------------- | ------------ |
| 0   | $x_0$ | $f_0$ | $f_1-f_0$  | $\Delta f_1 - \Delta f_0$ | ...          |
| 1   | $x_1$ | $f_1$ | $f_2-f_1$  | ...                       |              |
| 2   | $x_2$ | $f_2$ | ...        |                           |              |

2. Calcular $s = \dfrac{x - x_0}{h}$.
3. Evaluar la fórmula usando la **primera fila** de la tabla: $$P(x) = f_0 + s\,\Delta f_0 + \frac{s(s-1)}{2!}\Delta^2 f_0 + \cdots$$
**Ejemplo**: $f(x) = \ln(x)$, puntos en $x = 1.0, 1.1, 1.2, 1.3, 1.4$ ($h = 0.1$), evaluar en $x = 1.15$.
$s = \frac{1.15 - 1.0}{0.1} = 1.5$
$P(1.15) = 0 + 1.5(0.09531) + \frac{1.5(0.5)}{2}(-0.00830) + \cdots \approx 0.13976$

---
## <font color="#00b0f0">3. Newton hacia atrás</font>
**Cuándo usarlo**: datos **equiespaciados**, punto $x$ cerca del **final** de la tabla ($x_n$).

**Pasos**:
1. Construir tabla de diferencias $\nabla$ (mismo cálculo que $\Delta$, pero se llena de abajo hacia arriba): $$\nabla f_i = f_i - f_{i-1} \qquad \nabla^2 f_i = \nabla f_i - \nabla f_{i-1}$$
2. Calcular $s = \dfrac{x - x_n}{h}$ — nota: $s$ será **negativo** si $x < x_n$.
3. Usar la **última fila** de la tabla: $$P(x) = f_n + s\,\nabla f_n + \frac{s(s+1)}{2!}\nabla^2 f_n + \frac{s(s+1)(s+2)}{3!}\nabla^3 f_n + \cdots$$
**Diferencia clave con hacia adelante**: los signos en los factores cambian — $(s+1), (s+2), \ldots$ en vez de $(s-1), (s-2), \ldots$

---
## <font color="#00b0f0">4. Newton diferencias divididas</font>
**Cuándo usarlo**: datos con **cualquier espaciado** (la versión general de Newton).

**Pasos**:
1. Construir la tabla triangular columna por columna:

| $i$ | $x_i$ | Orden 0 | Orden 1 | Orden 2 |
|-----|--------|---------|---------|---------|
| 0   | $x_0$  | $f_0$   | $f[x_0,x_1]$ | $f[x_0,x_1,x_2]$ |
| 1   | $x_1$  | $f_1$   | $f[x_1,x_2]$ | |
| 2   | $x_2$  | $f_2$   | | |
$$f[x_i, x_{i+1}] = \frac{f_{i+1} - f_i}{x_{i+1} - x_i}$$$$f[x_i, x_{i+1}, x_{i+2}] = \frac{f[x_{i+1},x_{i+2}] - f[x_i,x_{i+1}]}{x_{i+2} - x_i}$$
2. Los coeficientes del polinomio son la **diagonal superior** (primera fila): $f[x_0],\ f[x_0,x_1],\ f[x_0,x_1,x_2], \ldots$

3. Evaluar acumulando el producto: $$P(x) = c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1) + c_3(x-x_0)(x-x_1)(x-x_2) + \cdots$$
**Truco de mano**: en el paso 3, calcular de adentro hacia afuera (Horner): $$P(x) = c_0 + (x-x_0)\bigl[c_1 + (x-x_1)\bigl[c_2 + (x-x_2)[\cdots]\bigr]\bigr]$$
---
## <font color="#00b0f0">5. Hermite</font>
**Cuándo usarlo**: conoces tanto $f(x_k)$ como $f'(x_k)$ en cada nodo. Da un polinomio de grado $2n+1$ que reproduce valores **y** derivadas.

**Pasos**:
1. Crear un vector $z$ duplicando cada nodo: $z_0 = z_1 = x_0,\ z_2 = z_3 = x_1, \ldots$
2. Inicializar la tabla de diferencias divididas con $f[z_i] = f(x_k)$.
3. Para cada par de nodos iguales, la diferencia dividida de orden 1 se define directamente como la derivada: $$f[z_{2k}, z_{2k+1}] = f'(x_k)$$
4. El resto de la tabla se completa con la fórmula normal de diferencias divididas.
5. Evaluar el polinomio resultante con la fórmula de Newton sobre los nodos $z$.

**Ejemplo** con 2 puntos ($n=1$): $(x_0, f_0, f'_0)$ y $(x_1, f_1, f'_1)$ → nodos $z_0=z_1=x_0,\ z_2=z_3=x_1$.

Tabla:
- $f[z_0] = f_0,\ f[z_1]=f_0,\ f[z_2]=f_1,\ f[z_3]=f_1$
- $f[z_0,z_1] = f'_0$ (definido así), $f[z_1,z_2] = \frac{f_1-f_0}{x_1-x_0}$, $f[z_2,z_3]=f'_1$
- Continuar la tabla normalmente.

### Ejemplo
Ejemplo completo paso a paso con números simples.
#### Datos del ejemplo

| $x_k$     | $f(x_k)$ | $f'(x_k)$ |
| --------- | -------- | --------- |
| $x_0 = 0$ | $1$      | $1$       |
| $x_1 = 1$ | $2$      | $3$       |

Tenemos $n=1$ (2 puntos), así que el polinomio resultante es de grado $2n+1 = 3$.

#### Paso 1 — Duplicar los nodos
Crear el vector $z$ repitiendo cada $x_k$: $$z_0 = z_1 = 0 \qquad z_2 = z_3 = 1$$
#### Paso 2 — Construir la tabla de diferencias divididas
**Columna 0** — copiar los valores de $f$:

| $i$ | $z_i$ | Ord. 0 |
| --- | ----- | ------ |
| 0   | 0     | 1      |
| 1   | 0     | 1      |
| 2   | 1     | 2      |
| 3   | 1     | 2      |
**Columna 1** — para nodos iguales se usa la derivada directamente: $$f[z_0, z_1] = f'(x_0) = 1 \quad \leftarrow \text{nodos iguales, derivada}$$$$f[z_1, z_2] = \frac{f[z_2]-f[z_1]}{z_2-z_1} = \frac{2-1}{1-0} = 1$$$$f[z_2, z_3] = f'(x_1) = 3 \quad \leftarrow \text{nodos iguales, derivada}$$
**Columna 2** — fórmula normal: $$f[z_0,z_1,z_2] = \frac{f[z_1,z_2] - f[z_0,z_1]}{z_2-z_0} = \frac{1-1}{1-0} = 0$$$$f[z_1,z_2,z_3] = \frac{f[z_2,z_3] - f[z_1,z_2]}{z_3-z_1} = \frac{3-1}{1-0} = 2$$
**Columna 3**:
$$f[z_0,z_1,z_2,z_3] = \frac{f[z_1,z_2,z_3] - f[z_0,z_1,z_2]}{z_3-z_0} = \frac{2-0}{1-0} = 2$$

**Tabla completa:**

| $i$ | $z_i$ | Ord. 0 | Ord. 1 | Ord. 2 | Ord. 3 |
| --- | ----- | ------ | ------ | ------ | ------ |
| 0   | 0     | **1**  | **1**  | **0**  | **2**  |
| 1   | 0     | 1      | 1      | 2      |        |
| 2   | 1     | 2      | 3      |        |        |
| 3   | 1     | 2      |        |        |        |

---
#### Paso 3 — Leer los coeficientes (primera fila en negrita)
$$c_0 = 1 \qquad c_1 = 1 \qquad c_2 = 0 \qquad c_3 = 2$$

#### Paso 4 — Armar el polinomio (Newton con nodos $z$)
$$H(x) = c_0 + c_1(x-z_0) + c_2(x-z_0)(x-z_1) + c_3(x-z_0)(x-z_1)(x-z_2)$$
$$H(x) = 1 + 1\cdot(x-0) + 0\cdot(x)(x) + 2\cdot(x)(x)(x-1)$$
$$\boxed{H(x) = 1 + x + 2x^2(x-1) = 2x^3 - 2x^2 + x + 1}$$

---
#### Verificación

| Condición           | Cálculo     | ¿Cumple? |
| ------------------- | ----------- | -------- |
| $H(0) = 1$          | $1$         | ✓        |
| $H(1) = 2$          | $2-2+1+1=2$ | ✓        |
| $H'(x) = 6x^2-4x+1$ | —           | —        |
| $H'(0) = 1$         | $1$         | ✓        |
| $H'(1) = 3$         | $6-4+1=3$   | ✓        |

El polinomio pasa por los dos puntos **y** tiene las dos derivadas correctas. Por eso con los mismos nodos Hermite es más preciso que Lagrange — exprime más información de cada punto.

---
## <font color="#00b0f0">6. Chebyshev</font>
**Cuándo usarlo**: cuando **tú eliges** los nodos y quieres minimizar el error máximo en $[a, b]$.

**Idea**: los nodos equiespaciados causan el fenómeno de Runge (error grande en extremos). Los nodos de Chebyshev son más densos en los extremos del intervalo, lo que iguala el error en todo el intervalo.

**Pasos**:
1. Calcular los nodos en $[a, b]$:
$$x_k = \frac{a+b}{2} + \frac{b-a}{2}\cos\!\left(\frac{2k+1}{2(n+1)}\,\pi\right), \quad k = 0, 1, \ldots, n$$
2. Evaluar $f$ en esos nodos.
3. Aplicar cualquier método de interpolación (Lagrange o Newton) sobre esos nodos.

**Ejemplo**: en $[-1, 1]$ con $n = 2$ (3 nodos):
$$x_0 = \cos\!\tfrac{\pi}{6} \approx 0.866 \qquad x_1 = \cos\!\tfrac{3\pi}{6} = 0 \qquad x_2 = \cos\!\tfrac{5\pi}{6} \approx -0.866$$

**Por qué funcionan**: el producto $\prod(x-x_i)$ con estos nodos tiene la cota más pequeña posible entre todos los conjuntos de $n+1$ puntos.

---
## <font color="#00b0f0">7. Splines cúbicos naturales</font>
**Cuándo usarlo**: muchos puntos y quieres una curva **suave** sin el problema de Runge. En vez de un polinomio de grado $n$, usa polinomios cúbicos en cada tramo.

**Pasos**:
1. Calcular los $h_i = x_{i+1} - x_i$.
2. Armar el sistema tridiagonal $(n-1)$ ecuaciones para $M_1, M_2, \ldots, M_{n-1}$: $$h_{i-1}M_{i-1} + 2(h_{i-1}+h_i)M_i + h_iM_{i+1} = 6\!\left(\frac{f_{i+1}-f_i}{h_i} - \frac{f_i-f_{i-1}}{h_{i-1}}\right)$$ con $M_0 = M_n = 0$ (condición natural).
3. Resolver el sistema (método de eliminación gaussiana o Thomas).
4. Calcular los coeficientes de cada tramo: $$a_i = f_i \qquad c_i = \frac{M_i}{2} \qquad d_i = \frac{M_{i+1}-M_i}{6h_i} \qquad b_i = \frac{f_{i+1}-f_i}{h_i} - \frac{h_i(2M_i+M_{i+1})}{6}$$
5. Para evaluar en $x$: identificar el intervalo $[x_i, x_{i+1}]$ que contiene $x$ y calcular $S_i(x)$.

### Splines Cúbicos
#### Cómo funciona el spline cúbico natural
**El problema con Lagrange y Newton**: un solo polinomio para muchos puntos → oscila fuera de control (Runge).

**La solución del spline**: en vez de un polinomio de grado $n$, usa $n$ cúbicos pegados: $$S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3 \quad \text{en cada } [x_i, x_{i+1}]$$

Para que los tramos "se vean suaves" al unirse, se impone:
- $S_i(x_{i+1}) = S_{i+1}(x_{i+1})$ — continua
- $S'_i(x_{i+1}) = S'_{i+1}(x_{i+1})$ — primera derivada continua
- $S''_i(x_{i+1}) = S''_{i+1}(x_{i+1})$ — segunda derivada continua ("curvatura" continua)

**Frontera natural**: $S''(x_0) = 0$ y $S''(x_n) = 0$ — la curva no tiene curvatura forzada en los extremos.

---
#### El truco clave: llamar $M_i = S''(x_i)$
Si defines $M_i$ como la segunda derivada en cada nodo, los 4 coeficientes se expresan directamente: $$a_i = f_i \qquad c_i = \frac{M_i}{2} \qquad d_i = \frac{M_{i+1}-M_i}{6h_i}$$$$b_i = \frac{f_{i+1}-f_i}{h_i} - \frac{h_i(2M_i+M_{i+1})}{6}$$
El único trabajo real es encontrar los $M_i$. Imponiendo continuidad de $S''$ aparece el **sistema tridiagonal** que ves en el formulario. Eso es lo que resuelve el **algoritmo de Thomas** — que el código implementa en `thomas_solve()`.

Una vez tienes los $M_i$, calculas los 4 coeficientes de cada tramo y evalúas en cualquier punto buscando en qué intervalo cae $x$.

---
#### Estructura del código
```
thomas_solve()   ← resuelve el sistema tridiagonal (O(n), muy rápido)
spline_fit()     ← construye el spline: calcula M_i y luego a,b,c,d
spline_eval()    ← evalúa S(x) en cualquier punto
Demo             ← mismos 6 nodos del quiz + comparación con Newton
```

Corre el archivo y al final verás que nuestra implementación coincide con `splinefun()` de R a 8 cifras decimales.

---
## <font color="#00b0f0">Comparación rápida</font>

| Método | Espaciado | Necesita | Ventaja |
|--------|-----------|----------|---------|
| Lagrange | Cualquiera | $f(x_k)$ | Fórmula directa, fácil de entender |
| Newton adelante | Uniforme | $f(x_k)$ | Evaluar cerca de $x_0$ |
| Newton atrás | Uniforme | $f(x_k)$ | Evaluar cerca de $x_n$ |
| Dif. divididas | Cualquiera | $f(x_k)$ | Agregar puntos sin recalcular |
| Hermite | Cualquiera | $f(x_k)$ y $f'(x_k)$ | Mayor precisión con mismos nodos |
| Chebyshev | Tú lo eliges | $f(x_k)$ | Mínimo error máximo posible |
| Splines | Cualquiera | $f(x_k)$ | Suavidad, sin Runge, muchos puntos |

**Todos producen el mismo polinomio único** (excepto Hermite y Splines que son constructions distintas).
