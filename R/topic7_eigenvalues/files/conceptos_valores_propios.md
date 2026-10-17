# Conceptos de apoyo — Valores Propios (Capítulo 7)

> Complemento de `guia_valores_propios.md`. Aquí no hay programas: solo las ideas
> que hay que tener claras para que la exposición no se sienta como una receta.
> Todos los ejemplos numéricos de este archivo están verificados en R.

**Índice**

1. [Por qué el problema es $A\mathbf{x} = \lambda\mathbf{x}$ (repaso de álgebra lineal)](#1-por-qué-el-problema-es-axlambdax)
2. [Qué es $\omega$ y por qué aparece en la motivación](#2-qué-es-ω-y-por-qué-aparece-la-vibración)
3. [Similaridad, Householder y QR](#3-similaridad-householder-y-qr)
4. [Matriz no simétrica: ¿valores propios complejos?](#4-matriz-no-simétrica-valores-propios-complejos)
5. [Cómo $N+1$ puntos "determinan" un polinomio (y qué NO hacen)](#5-cómo-n1-puntos-determinan-un-polinomio)
6. [Diferencias hacia adelante de Newton y coeficientes de Markov](#6-diferencias-de-newton-y-coeficientes-de-markov)
7. [Qué es una matriz tridiagonal](#7-qué-es-una-matriz-tridiagonal)
8. [Cómo buscar videos en YouTube](#8-cómo-buscar-videos-en-youtube)

---

# 1. Por qué el problema es $A\mathbf{x} = \lambda\mathbf{x}$

## La idea en una frase

> Una matriz $A$ **transforma** vectores: les cambia el tamaño *y* la dirección.
> Casi siempre cambia ambas cosas. Un **vector propio** es una dirección especial
> que $A$ **no gira**: solo la estira o la encoge. El factor de estiramiento es el
> **valor propio** $\lambda$.

$$A\mathbf{x} = \lambda\mathbf{x}$$

"Aplicar $A$ a $\mathbf{x}$ da lo mismo que multiplicar $\mathbf{x}$ por un número."

## Diferencia con el sistema de siempre

| | Sistema lineal | Problema de valores propios |
|---|---|---|
| Ecuación | $A\mathbf{x} = \mathbf{y}$ | $A\mathbf{x} = \lambda\mathbf{x}$ |
| Datos | $A$ y $\mathbf{y}$ | solo $A$ |
| Incógnitas | $\mathbf{x}$ | $\lambda$ **y** $\mathbf{x}$ |
| Lado derecho | un vector conocido | $\lambda\mathbf{x}$, desconocido |
| Cuántas soluciones | una (si $\det A \ne 0$) | $N$ valores de $\lambda$, cada uno con su dirección |

## De dónde sale $\det(A - \lambda I) = 0$

Pasamos todo a un lado:

$$A\mathbf{x} - \lambda\mathbf{x} = \mathbf{0} \;\Longrightarrow\; (A - \lambda I)\,\mathbf{x} = \mathbf{0}$$

(Se escribe $\lambda I$ y no $\lambda$ a secas porque no se puede restar un número a una
matriz; $\lambda I$ es la matriz que tiene $\lambda$ en la diagonal.)

Eso es un **sistema homogéneo**. Siempre tiene la solución trivial $\mathbf{x} = \mathbf{0}$,
pero esa no sirve (un vector propio es, por definición, no nulo). Y aquí viene el
repaso:

| Si $\det(A - \lambda I)$ … | entonces … |
|---|---|
| $\ne 0$ | la matriz es invertible, la **única** solución es $\mathbf{x}=\mathbf{0}$ ✗ |
| $= 0$ | la matriz es singular, **hay infinitas** soluciones no triviales ✓ |

Por eso buscamos los $\lambda$ que hacen **singular** a $A - \lambda I$:

$$\boxed{\det(A - \lambda I) = 0}$$

Es justo lo que pasa en el ejemplo del inicio del tema:

$$3x - 2y + z = 0,\quad x + y + 2z = 0,\quad 4x - y + 3z = 0$$

Su determinante es $0$ (la tercera ecuación es la suma de las dos primeras), y por eso
tiene soluciones no triviales: $(x,y,z) = \beta(1, 1, -1)$ para cualquier $\beta$.

> **Fíjate en el $\beta$:** el vector propio no es único, cualquier múltiplo también sirve.
> Si $A\mathbf{x} = \lambda\mathbf{x}$, entonces $A(5\mathbf{x}) = \lambda(5\mathbf{x})$.
> Lo que importa es la **dirección**, no el tamaño.

## Ejemplo completo $2\times2$

$$A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$$

**Paso 1 — polinomio característico:**

$$\det(A - \lambda I) = (2-\lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3)$$

Valores propios: $\lambda_1 = 3$, $\lambda_2 = 1$.

**Paso 2 — un vector propio por cada $\lambda$:**

| $\lambda$ | $(A - \lambda I)\mathbf{x} = \mathbf{0}$ | Vector propio |
|:---:|---|:---:|
| $3$ | $\begin{bmatrix} -1 & 1 \\ 1 & -1 \end{bmatrix}\mathbf{x} = \mathbf{0} \Rightarrow y = x$ | $(1,\ 1)$ |
| $1$ | $\begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix}\mathbf{x} = \mathbf{0} \Rightarrow y = -x$ | $(1,\ -1)$ |

**Paso 3 — comprobar:**

$$A\begin{bmatrix}1\\1\end{bmatrix} = \begin{bmatrix}3\\3\end{bmatrix} = 3\begin{bmatrix}1\\1\end{bmatrix} \checkmark
\qquad
A\begin{bmatrix}1\\-1\end{bmatrix} = \begin{bmatrix}1\\-1\end{bmatrix} = 1\begin{bmatrix}1\\-1\end{bmatrix} \checkmark$$

**Lo que hace $A$ geométricamente:** estira ×3 en la dirección de la diagonal $(1,1)$ y
deja quieta la dirección $(1,-1)$. Cualquier otro vector sí cambia de dirección:
$(1,0) \mapsto (2,1)$.

## Para refrescar: diagonalización

Si juntamos los vectores propios como columnas de $P$ y los valores propios en una
diagonal $D$:

$$A = P\,D\,P^{-1}, \qquad P = \begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}, \quad D = \begin{bmatrix} 3 & 0 \\ 0 & 1 \end{bmatrix}$$

(se puede verificar multiplicando: $PDP^{-1}$ devuelve $A$). Es el $P$ y el $P^{-1}$ del
curso de Cálculo; en la sección 3 se conecta con la lámina.

## Chequeos que siempre conviene hacer

$$\lambda_1 + \lambda_2 + \cdots + \lambda_N = \operatorname{traza}(A), \qquad \lambda_1\lambda_2\cdots\lambda_N = \det(A)$$

En el ejemplo: $3 + 1 = 4 = 2 + 2$ ✓ y $3 \cdot 1 = 3 = 4 - 1$ ✓.

## Una propiedad que se usa mucho: simétrica ⟹ todo es bonito

Si $A$ es **simétrica** ($A = A^T$) y real:
- todos los valores propios son **reales**,
- los vectores propios de valores distintos son **ortogonales**
  (en el ejemplo, $(1,1)\cdot(1,-1) = 0$).

Esa es la razón de que el Programa 7-2 solo funcione con simétricas.

---

# 2. Qué es $\omega$ y por qué aparece la vibración

## Respuesta corta

**Sí, $\omega$ es una velocidad angular** — más precisamente, la **frecuencia angular**
de la oscilación, en rad/s. Aquí no hay nada girando; es la "velocidad" con la que
avanza el ángulo dentro de $\cos(\omega t)$.

| Magnitud | Símbolo | Relación | Unidad |
|---|:---:|---|:---:|
| Frecuencia angular | $\omega$ | — | rad/s |
| Frecuencia | $f$ | $f = \omega/2\pi$ | Hz (ciclos/s) |
| Período | $T$ | $T = 2\pi/\omega$ | s |

> **Intuición:** imagina un punto girando en un círculo con velocidad angular $\omega$.
> Su sombra proyectada sobre una pared sube y baja como $\cos(\omega t)$. Una masa en un
> resorte se mueve exactamente así.

## Por qué aparece $e^{i\omega t}$

El sistema de masas y resortes obedece (para una masa suelta con resorte $k$):

$$m\,y''(t) = -k\,y(t)$$

Es una ecuación diferencial: relaciona $y$ con su derivada segunda. El truco es
**proponer** una solución que oscile a una sola frecuencia:

$$y(t) = e^{i\omega t}f, \qquad e^{i\omega t} = \cos\omega t + i\sin\omega t$$

¿Por qué esa forma? Porque derivarla dos veces da el mismo vector multiplicado por una
constante:

$$y'' = (i\omega)^2\,e^{i\omega t}f = -\omega^2\,y$$

Sustituyendo, **desaparece el tiempo** y la ecuación diferencial se vuelve algebraica:

$$m(-\omega^2)f = -kf \;\Longrightarrow\; \frac{k}{m}f = \omega^2 f$$

Con varias masas, $k/m$ pasa a ser una matriz $A = K/m$, y queda:

$$\boxed{A\mathbf{f} = \omega^2\,\mathbf{f}} \qquad\Longleftrightarrow\qquad A\mathbf{f} = \lambda\mathbf{f}, \quad \lambda = \omega^2$$

> **Ahí está la conexión.** Pedir que *todas* las masas oscilen con la misma
> frecuencia $\omega$ es pedir que $\mathbf{f}$ sea un vector propio. El movimiento real
> es la parte real: $y_j(t) = f_j\cos(\omega t)$.

## Qué significa cada cosa

- **Valor propio $\lambda = \omega^2$** → cuadrado de una **frecuencia natural** del sistema.
- **Vector propio $\mathbf{f}$** → **forma de vibrar** (el "modo"): cuánto se mueve cada
  masa respecto a las otras.

## Con números ($k = m = 1$, tres masas, cuatro resortes)

$$A = \begin{bmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{bmatrix}$$

| Modo | $\lambda = \omega^2$ | $\omega$ (rad/s) | $f$ (Hz) | $T$ (s) | Forma |
|:---:|:---:|:---:|:---:|:---:|---|
| 1 | $2-\sqrt2 = 0.586$ | 0.765 | 0.122 | 8.21 | las tres masas en fase |
| 2 | $2$ | 1.414 | 0.225 | 4.44 | la del centro quieta, las otras en oposición |
| 3 | $2+\sqrt2 = 3.414$ | 1.848 | 0.294 | 3.40 | alternadas (arriba–abajo–arriba) |

## Por qué importa que $\lambda$ sea positivo

El tema dice que estos problemas sirven para "estudiar la **estabilidad dinámica**".
Esa es la razón:

| $\lambda = \omega^2$ | $\omega$ | $e^{i\omega t}$ | Comportamiento |
|---|---|---|---|
| $> 0$ | real | oscila | **estable** (vibra) |
| $< 0$ | imaginario | crece o decae exponencialmente | **inestable** |

Y como $K$ es **simétrica** (cada resorte tira igual en ambos sentidos), los $\lambda$
salen **reales**, y por eso no hay sorpresas con números complejos.

> **Detalle de signos:** en la matriz, la diagonal es positiva y fuera de la diagonal va
> $-k/m$. Sale de pasar al otro lado de la ecuación:
> $-\omega^2f_1 = -\tfrac1m(k_{01}+k_{12})f_1 + \tfrac1m k_{12}f_2$ → multiplicando por
> $-1$, el término del vecino queda $-k_{12}/m$.

---

# 3. Similaridad, Householder y QR

## La pregunta: ¿son los $P$ y $P^{-1}$ del curso de Cálculo?

**Sí, es la misma idea** (cambio de base: $P A P^{-1}$ o $P^{-1}AP$). Lo que cambia es
**de dónde sale $P$** y **qué se busca**:

| | Diagonalización / Jordan (Cálculo) | Householder y QR (este tema) |
|---|---|---|
| Forma | $A = PDP^{-1}$ | $B = PAP^{-1}$ |
| ¿De dónde sale $P$? | de los **vectores propios** | de las **entradas de $A$**, sin conocerlos |
| Resultado | diagonal $D$ en un paso | una forma **más simple** (tridiagonal / triangular) |
| Problema | para construir $P$ hay que *haber resuelto* el problema | **no** requiere conocer los vectores propios |

> **La idea de fondo:** diagonalizar de golpe es lo ideal, pero requiere conocer lo que
> buscamos. Así que se hace por **etapas**: transformaciones sencillas que no cambian los
> valores propios y que dejan la matriz cada vez más simple, hasta que los valores
> propios se pueden leer o calcular fácilmente.

## Qué es una transformación de similaridad

Dos matrices $A$ y $B$ son **semejantes** si existe $P$ invertible con

$$B = P\,A\,P^{-1}$$

**Teorema:** semejantes ⟹ **mismos valores propios**.

**Por qué** (3 líneas):

$$\det(B - \lambda I) = \det\bigl(PAP^{-1} - \lambda PP^{-1}\bigr) = \det\bigl(P(A-\lambda I)P^{-1}\bigr)
= \underbrace{\det P\cdot\det P^{-1}}_{=1}\cdot\det(A - \lambda I)$$

Mismo polinomio característico ⟹ mismas raíces.

**Los vectores propios sí cambian**: si $A\mathbf{x} = \lambda\mathbf{x}$, entonces
$B(P\mathbf{x}) = \lambda(P\mathbf{x})$. El vector propio de $B$ es $P\mathbf{x}$.

## Ejemplo con números

$$A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}, \qquad P = \begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}, \qquad P^{-1} = \begin{bmatrix} 1 & -1 \\ 0 & 1 \end{bmatrix}$$

$$B = PAP^{-1} = \begin{bmatrix} 3 & 0 \\ 1 & 1 \end{bmatrix}$$

| | $A$ | $B$ |
|---|:---:|:---:|
| Valores propios | $3,\ 1$ | $3,\ 1$ ✓ |
| Traza | 4 | 4 ✓ |
| Determinante | 3 | 3 ✓ |
| Vector propio de $\lambda=3$ | $(1,1)$ | $P(1,1) = (2,1)$ ✓ |

$B$ es triangular, así que sus valores propios **se leen en la diagonal**. Ese es el
objetivo de todos los métodos de este capítulo.

## Householder

**Qué es:** una matriz $P$ que es un **espejo** (reflexión):

$$P = I - \frac{\mathbf{u}\mathbf{u}^T}{h}, \qquad h = \tfrac12\mathbf{u}^T\mathbf{u}$$

**Propiedades** (las del tema): $P = P^T = P^{-1}$, es decir, $PP = I$. Reflejar dos
veces devuelve el vector original, y una reflexión **no cambia longitudes**.

**Para qué se usa:** se elige $\mathbf{u}$ a medida para que el espejo **convierta en
ceros** los elementos de una columna que están debajo de la subdiagonal. Repetido en
cada columna ($N-2$ veces), una matriz simétrica queda **tridiagonal**.

Cada paso es $A \to PAP$ (con $P^{-1} = P$), una similaridad: **los valores propios no
cambian**. En el ejemplo del tema:

| | $A$ original | tridiagonal $T$ |
|---|:---:|:---:|
| Traza | $3 + 7 + 0 = 10$ | $3 + 1.3529 + 5.6471 = 10$ ✓ |

## QR

**Qué es:** factorizar $A = QR$ con $Q$ **ortogonal** ($Q^{-1} = Q^T$) y $R$ **triangular
superior**. Después se multiplica en orden inverso:

$$A_{k+1} = R_kQ_k = Q_k^{-1}\,(Q_kR_k)\,Q_k = Q_k^{-1}A_kQ_k$$

La segunda igualdad muestra que cada paso es una **similaridad** (los valores propios
no cambian). Y repetido muchas veces, $A_k$ **converge a una matriz triangular superior**
(por bloques, si hay complejos), con los valores propios en la diagonal.

> **Resumen de los tres:**
> - **Similaridad** = la regla que garantiza que no se pierde nada.
> - **Householder** = una forma concreta de $P$ para *tridiagonalizar* en un número finito de pasos.
> - **QR** = una forma concreta de iterar que *converge* a triangular.

---

# 4. Matriz no simétrica: ¿valores propios complejos?

## Respuesta

| Matriz real… | Valores propios |
|---|---|
| **simétrica** | **siempre reales** (garantizado) |
| **no simétrica** | *pueden* ser reales **o** complejos |

"Puede", no "tiene". Y cuando aparecen complejos, siempre salen en **pares conjugados**
$a \pm bi$ (porque el polinomio característico tiene coeficientes reales).

## Ejemplos

| Matriz | No simétrica | $\lambda$ | ¿Complejos? |
|---|:---:|---|:---:|
| $\begin{bmatrix} 1 & 2 \\ 0 & 3 \end{bmatrix}$ | ✓ | $1,\ 3$ | no |
| $\begin{bmatrix} 0 & -1 \\ 1 & 0 \end{bmatrix}$ | ✓ | $\pm i$ | **sí** |

**La segunda es una rotación de $90^\circ$:**
$\det(A - \lambda I) = \lambda^2 + 1 = 0 \Rightarrow \lambda = \pm i$.

> **Geometría:** una rotación gira *todas* las direcciones, así que no hay ningún vector
> real que conserve su dirección. Por eso no hay valores propios reales: los complejos
> son la manera del álgebra de decir "esto gira".

## Dos hechos útiles

- Polinomio de grado **impar** con coeficientes reales ⟹ al menos **una raíz real**.
  Entonces una matriz real de orden impar siempre tiene algún valor propio real.
- En el tema: el polinomio de grado 5 con la matriz de Frobenius (**no simétrica**) dio
  $-0.9086$ y los pares $0.1404 \pm 2.7314i$ y $0.4139 \pm 0.3506i$. Un real y dos pares
  complejos: $1 + 2 + 2 = 5$.

## Por qué cada programa necesita lo que necesita

| Programa | Matriz | Por qué |
|---|---|---|
| **7-1** (interpolación) | cualquiera | solo necesita calcular determinantes |
| **7-2** (Householder + bisección) | **simétrica** | Householder da tridiagonal solo si hay simetría; la bisección busca raíces **reales** en intervalos |
| QR | cualquiera | por eso devuelve bloques $2\times2$ para los complejos |

---

# 5. Cómo $N+1$ puntos "determinan" un polinomio

## La confusión más común (y es muy razonable)

> **Evaluar el polinomio en $N+1$ puntos NO encuentra sus raíces.**
> Lo que hace es **determinar el polinomio completo**: sus $N+1$ coeficientes.

Son dos problemas distintos:

| Problema | Entrada | Salida |
|---|---|---|
| **Determinar** el polinomio (Programa 7-1) | $N+1$ valores $f(\lambda_0), \dots, f(\lambda_N)$ | los coeficientes $c_0, \dots, c_N$ |
| **Hallar las raíces** (paso siguiente) | los coeficientes | los $\lambda$ con $f(\lambda) = 0$ |

## La analogía de siempre

| Grado | Forma | Coeficientes | Puntos que lo fijan |
|:---:|---|:---:|:---:|
| 1 | recta $c_0 + c_1\lambda$ | 2 | 2 |
| 2 | parábola $c_0 + c_1\lambda + c_2\lambda^2$ | 3 | 3 |
| $N$ | $c_0 + \cdots + c_N\lambda^N$ | $N+1$ | $N+1$ |

Por dos puntos pasa **una sola** recta; por tres, **una sola** parábola.

## Por qué: es un sistema lineal

Un polinomio de grado $N$ tiene $N+1$ coeficientes desconocidos. Cada punto da una
ecuación:

$$c_0 + c_1\lambda_j + c_2\lambda_j^2 + \cdots + c_N\lambda_j^N = f(\lambda_j), \qquad j = 0, \dots, N$$

Son $N+1$ ecuaciones con $N+1$ incógnitas. La matriz es:

$$\begin{bmatrix} 1 & \lambda_0 & \lambda_0^2 & \cdots \\ 1 & \lambda_1 & \lambda_1^2 & \cdots \\ \vdots & & & \\ 1 & \lambda_N & \lambda_N^2 & \cdots \end{bmatrix}$$

> **¡Es una matriz de Vandermonde!** La misma del Programa 5-1, transpuesta. Es invertible
> si los $\lambda_j$ son distintos. Por eso la solución existe y es única.

**Otra forma de verlo (unicidad):** si dos polinomios de grado $\le N$ coincidieran en
$N+1$ puntos, su diferencia tendría grado $\le N$ y $N+1$ raíces. Eso solo lo cumple el
polinomio cero. Luego son el mismo.

## Qué alimenta el programa y qué determina

| | Qué es | Quién lo elige |
|---|---|---|
| **Puntos** $\lambda_j = j\,\Delta\lambda$ | dónde evaluamos | nosotros (equiespaciados, empezando en 0) |
| **Valores** $f(\lambda_j)$ | $\det(A - \lambda_jI)$ | se **calculan** (Gauss, cap. 6) |
| **Incógnitas** | los $N+1$ coeficientes $c_i$ | las halla el método |

Entonces sí: se alimenta con las **imágenes** $f(\lambda_j)$, que son los determinantes.
Los puntos $\lambda_j$ son arbitrarios (no son valores propios, no los conocemos); lo
único que hace falta es que sean distintos.

> **Ojo:** en el ejemplo del tema, $f(1) = 0$. Eso es una *casualidad* (resulta que
> $\lambda = 1$ es valor propio). En general, los valores $f(\lambda_j)$ son números
> cualesquiera, distintos de cero.

## Ejemplo mínimo: $A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$ (grado 2 → 3 puntos)

**Paso 1 — evaluar** en $\lambda = 0, 1, 2$ (con $\Delta\lambda = 1$):

| $\lambda$ | $A - \lambda I$ | $f(\lambda) = \det$ |
|:---:|:---:|:---:|
| 0 | $\begin{bmatrix} 2&1\\1&2 \end{bmatrix}$ | $3$ |
| 1 | $\begin{bmatrix} 1&1\\1&1 \end{bmatrix}$ | $0$ |
| 2 | $\begin{bmatrix} 0&1\\1&0 \end{bmatrix}$ | $-1$ |

**Paso 2 — reconstruir** (con las diferencias de Newton, sección 6):

$$f(s) = 3 - 3s + 2\binom{s}{2} = 3 - 3s + (s^2 - s) = s^2 - 4s + 3$$

**Paso 3 — recién ahora, las raíces:** $\lambda^2 - 4\lambda + 3 = 0 \Rightarrow \lambda = 1, 3$.

Nunca evaluamos en $\lambda = 3$; el valor propio $3$ apareció **después**, al resolver el
polinomio.

## El flujo completo

```
matriz A
   │  N+1 determinantes (Gauss)          ← Programa 7-1, parte 1
   ▼
valores f(λ₀), f(λ₁), ..., f(λ_N)
   │  diferencias de Newton + Markov     ← Programa 7-1, parte 2
   ▼
coeficientes c₀, c₁, ..., c_N            ← resultado del 7-1
   │  Frobenius + QR  (o bisección)      ← paso aparte
   ▼
raíces = valores propios
```

---

# 6. Diferencias de Newton y coeficientes de Markov

## Repaso: la fórmula de Newton hacia adelante

Si tenemos valores $f_0, f_1, \dots, f_N$ en puntos **equiespaciados**
$\lambda_j = \lambda_0 + j\,\Delta\lambda$, definimos $s = (\lambda - \lambda_0)/\Delta\lambda$ y:

$$f(s) = f_0 + \Delta f_0\,s + \Delta^2 f_0\,\frac{s(s-1)}{2!} + \Delta^3 f_0\,\frac{s(s-1)(s-2)}{3!} + \cdots$$

o, en forma compacta:

$$f(s) = \sum_{k=0}^{N}\Delta^k f_0\,\binom{s}{k}, \qquad \binom{s}{k} = \frac{s(s-1)\cdots(s-k+1)}{k!}$$

donde $\Delta f_j = f_{j+1} - f_j$ y $\Delta^k$ repite la operación $k$ veces.

**La tabla** (es la del tema 2): cada columna es la resta de la anterior.

| $s$ | $f$ | $\Delta f$ | $\Delta^2 f$ |
|:---:|:---:|:---:|:---:|
| 0 | **3** | **−3** | **2** |
| 1 | 0 | −1 | |
| 2 | −1 | | |

(ejemplo $2\times2$ de la sección anterior). **Solo se usa la primera fila**: $3,\ -3,\ 2$.

## Tu duda: ¿los $\Delta^k f_0$ son los coeficientes del polinomio?

**Casi.** Son los coeficientes del polinomio, pero **en otra base**.

Un mismo polinomio se puede escribir de dos maneras:

| Base | Cómo se ve | Coeficientes |
|---|---|---|
| **Potencias** (la que queremos) | $c_0 + c_1 s + c_2 s^2 + \cdots$ | $c_i$ |
| **Newton** (la que sale de la tabla) | $\Delta^0f_0 + \Delta^1f_0\,s + \Delta^2f_0\binom{s}{2} + \cdots$ | $\Delta^kf_0$ |

Son **el mismo polinomio** escrito con otros "ladrillos". La tabla de diferencias nos
da los coeficientes de Newton fácilmente (solo restas), pero lo que se necesita es el
polinomio en potencias de $\lambda$ para después buscar sus raíces.

## Entonces, ¿qué son los coeficientes de Markov?

> **Son el diccionario entre las dos bases:** los coeficientes de cada $\binom{s}{k}$
> escrito como polinomio en potencias de $s$.

| $k$ | $\binom{s}{k}$ | Desarrollado | Coeficientes de Markov (de $s^1$ en adelante) |
|:---:|---|---|---|
| 1 | $s$ | $s$ | $1$ |
| 2 | $\dfrac{s(s-1)}{2}$ | $\dfrac{s^2 - s}{2}$ | $-0.5,\ \ 0.5$ |
| 3 | $\dfrac{s(s-1)(s-2)}{6}$ | $\dfrac{s^3 - 3s^2 + 2s}{6}$ | $0.3333,\ -0.5,\ 0.1667$ |
| 4 | $\dfrac{s(s-1)(s-2)(s-3)}{24}$ | $\dfrac{s^4 - 6s^3 + 11s^2 - 6s}{24}$ | $-0.25,\ 0.4583,\ -0.25,\ 0.0417$ |

Son **exactamente** las filas que imprime el Programa 7-1 bajo `COEFICIENTES DE MARKOV`:

```
 -0.5000000    0.5000000
  0.3333333   -0.5000000    0.1666667
 -0.2500000    0.4583333   -0.2500000    0.0416667
```

> **Sobre el nombre:** así los llama el libro. En otros textos aparecen como los
> **números de Stirling de primera especie** divididos entre $k!$. Lo importante es qué
> son (los coeficientes de $\binom{s}{k}$ en potencias), no el nombre.

## Cómo se usan (el ejemplo $2\times2$ completo)

$f(s) = 3\cdot\mathbf{1} \;+\; (-3)\cdot s \;+\; 2\cdot\dfrac{s^2 - s}{2}$

Se suma **potencia por potencia**:

| Potencia | Aporta | Total |
|:---:|---|:---:|
| $s^0$ | $3$ | $3$ |
| $s^1$ | $-3$ (de $\Delta f_0$) $+\ 2(-0.5)$ (de $\Delta^2 f_0$) | $-4$ |
| $s^2$ | $2(0.5)$ (de $\Delta^2 f_0$) | $1$ |

$$f(s) = s^2 - 4s + 3$$

Como aquí $\Delta\lambda = 1$, $s = \lambda$ y ya está. Si $\Delta\lambda \ne 1$, hay un
último paso: $s^i = \lambda^i/\Delta\lambda^{\,i}$, es decir, **dividir $c_i$ entre
$\Delta\lambda^i$** (la línea `CC <- CC / DE^(0:N)` del programa).

## Cómo lo hace el programa (la línea que asusta)

```r
MV[L + 1] <- (MV[L] - (K - 1) * MV[L + 1]) / K
```

Es solo la relación

$$\binom{s}{K} = \binom{s}{K-1}\cdot\frac{s - (K-1)}{K}$$

aplicada a los coeficientes. Multiplicar por $s$ **corre** los coeficientes una posición
(`MV[L]`); multiplicar por $-(K-1)$ los **escala** (`(K-1) * MV[L+1]`); y se divide entre
$K$. Así se pasa de la fila $K-1$ a la fila $K$ sin recalcular desde cero.

## Una comprobación que da confianza

Evaluemos la fórmula de Newton en los propios puntos de la tabla ($s = 0, 1, 2$):

| $s$ | $f_0 + \Delta f_0\,s + \Delta^2f_0\binom{s}{2}$ | Resultado | ¿Coincide? |
|:---:|---|:---:|:---:|
| 0 | $3 + 0 + 0$ | $3$ | $f_0 = 3$ ✓ |
| 1 | $3 - 3 + 0$ | $0$ | $f_1 = 0$ ✓ |
| 2 | $3 - 6 + 2$ | $-1$ | $f_2 = -1$ ✓ |

La fórmula **pasa por todos los puntos**, y como el grado es $N$ con $N+1$ puntos, no
puede ser otro polinomio.

---

# 7. Qué es una matriz tridiagonal

## Definición

Una matriz **tridiagonal** solo tiene elementos no nulos en **tres diagonales**: la
principal, la de arriba y la de abajo. Todo lo demás es cero.

$$T = \begin{bmatrix}
a_1 & b_1 & 0 & 0 \\
b_1 & a_2 & b_2 & 0 \\
0 & b_2 & a_3 & b_3 \\
0 & 0 & b_3 & a_4
\end{bmatrix}$$

- $a_i$: diagonal principal ($N$ elementos).
- $b_i$: diagonal superior e inferior ($N-1$ elementos, **iguales** si es simétrica).

El del ejemplo del tema (el que sale de Householder):

$$\begin{bmatrix} 3 & -4.1231 & 0 \\ -4.1231 & 1.3529 & 3.4118 \\ 0 & 3.4118 & 5.6471 \end{bmatrix}$$

Los dos ceros de las esquinas son lo que Householder logró.

## Por qué se busca llegar a esa forma

| Ventaja | Detalle |
|---|---|
| **Cabe en poco espacio** | se guardan $3N - 2$ números y no $N^2$ |
| **Determinante barato** | recurrencia de 3 términos, $O(N)$ en vez de $O(N^3)$ |
| **Todas las submatrices gratis** | $p_i$ es el polinomio de la submatriz $i\times i$ |
| **Entrelazamiento** | permite encerrar cada raíz en un intervalo (bisección) |

La recurrencia de la que sale todo:

$$p_i(\lambda) = (a_i - \lambda)\,p_{i-1}(\lambda) - b_{i-1}^2\,p_{i-2}(\lambda)$$

(sale de desarrollar el determinante por la última fila: como solo hay dos elementos no
nulos, solo aparecen dos términos).

## Dónde aparecen en la naturaleza

Una matriz es tridiagonal cuando **cada elemento solo interactúa con sus vecinos
inmediatos**:

- **Masas en cadena:** cada masa solo siente a las dos vecinas (el ejemplo del tema).
- **Diferencias finitas** (capítulo 5): la fórmula de la segunda derivada,

$$f''_i \approx \frac{f_{i-1} - 2f_i + f_{i+1}}{h^2}$$

  produce una matriz con $-2$ en la diagonal y $1$ a los lados. Es la misma $K$ de las
  masas y resortes (salvo el signo): **la matriz de vibraciones es una segunda derivada
  discretizada**.
- **Splines cúbicos** (tema 2): el sistema que se resuelve es tridiagonal.

## Parientes

| Forma | Ceros | Cuándo sale |
|---|---|---|
| **Tridiagonal** | todo salvo 3 diagonales | simétrica + Householder |
| **Hessenberg superior** | todo debajo de la 1.ª subdiagonal | no simétrica + Householder |
| **Triangular superior** | todo debajo de la diagonal | tras la iteración QR |

> **Una tridiagonal simétrica es un caso particular de Hessenberg.** El texto lo dice:
> "la matriz tridiagonal es un caso especial de la matriz de Hessenberg".

---

# 8. Cómo buscar videos en YouTube

## Primero, una advertencia honesta

Desde aquí no puedo buscar en YouTube ni comprobar qué videos siguen disponibles.
Lo que sigue son **estrategias de búsqueda** y canales que conozco; verifica
cuáles te sirven.

## Truco 1 — buscar por el nombre del método, no por el del programa

El libro es *Applied Numerical Methods with Software* (Nakamura). Casi nadie hace videos
del "Programa 7-2". Lo que sí existe son videos de cada **método**. Traduce cada
programa a términos genéricos:

| Lo que dice tu material | Qué buscar (español) | Qué buscar (inglés) |
|---|---|---|
| Valores propios | `valores y vectores propios explicación` | `eigenvalues and eigenvectors` |
| Polinomio característico | `polinomio característico ejemplo resuelto` | `characteristic polynomial` |
| Método de las potencias | `método de las potencias valor propio dominante` | `power method eigenvalue` |
| Potencias inversas | `método de la potencia inversa` | `inverse power iteration` |
| Householder | `transformación de Householder tridiagonal` | `Householder tridiagonalization` |
| Bisección sobre tridiagonal | `bisección valores propios tridiagonal` | `Sturm sequence bisection eigenvalues` |
| Iteración QR | `algoritmo QR valores propios` | `QR algorithm eigenvalues` |
| Matriz de Frobenius / compañera | `matriz compañera raíces polinomio` | `companion matrix roots polynomial` |
| Diferencias de Newton | `interpolación de Newton diferencias progresivas` | `Newton forward difference interpolation` |
| Vibración con masas y resortes | `modos normales masas resortes valores propios` | `normal modes coupled oscillators eigenvalues` |

## Truco 2 — palabras que mejoran el resultado

- Añade **`ejemplo resuelto`**, **`paso a paso`** o **`worked example`**: filtra las
  clases solo teóricas.
- Para entender la **idea** (no la receta), busca **`intuición geométrica`** o
  **`visual`**.
- Para ver el **código**, añade `R`, `Python`, `MATLAB` o `código`. Si encuentras la
  versión en Python, es casi línea a línea igual a R.
- Usa un ejemplo con tu misma matriz $3\times3$ o $2\times2$: a veces conviene buscar
  directamente `[[2,1],[1,2]] eigenvalues`.

## Truco 3 — canales (de memoria, verifica)

| Canal | Idioma | Para qué sirve |
|---|---|---|
| **3Blue1Brown** — serie *Essence of Linear Algebra* | inglés (con subtítulos en español) | el mejor para la **intuición** de vectores propios y transformaciones. Mira el capítulo de "Eigenvectors and eigenvalues" |
| **MIT OpenCourseWare** — Gilbert Strang, 18.06 | inglés | lecciones completas de álgebra lineal, incluida la de valores propios |
| **Holistic Numerical Methods** (Autar Kaw) | inglés | métodos numéricos con ejemplo resuelto, muy parecido al formato de Nakamura |
| **Steve Brunton** | inglés | valores propios, SVD, visual y rápido |
| **Julioprofe** | español | álgebra lineal con ejemplos resueltos, estilo pizarrón |

## Truco 4 — un orden de visionado para este tema

1. **Intuición** ($A\mathbf{x} = \lambda\mathbf{x}$ y para qué sirve): la serie de 3Blue1Brown.
2. **Cálculo a mano** (polinomio característico, vectores propios): un video en español
   con ejemplo $2\times2$ y $3\times3$.
3. **Métodos numéricos**: potencias primero (es el más simple), luego QR; Householder
   al final porque es el más técnico.

## Truco 5 — cómo aprovechar los videos

- **Velocidad 1.25×–1.5×** en la parte teórica; **1×** cuando resuelven el ejemplo.
- Pausa antes de que terminen el ejemplo e intenta cerrarlo tú. Después compara.
- **Usa los videos para entender, no para memorizar**: lo que se evalúa en la exposición
  es poder explicar *por qué*, y eso lo tienes en `guia_valores_propios.md` y aquí.
- Si un video usa una notación distinta ($\lambda$ vs. $\omega$, $PAP^{-1}$ vs.
  $P^{-1}AP$), no te asustes: es el mismo concepto.

---

# Mini-glosario

| Término | Significa |
|---|---|
| **Valor propio** $\lambda$ | factor por el que $A$ estira un vector propio |
| **Vector propio** $\mathbf{x}$ | dirección que $A$ no gira |
| **Polinomio característico** | $\det(A - \lambda I)$, de grado $N$ |
| **Singular** | $\det = 0$, no invertible |
| **Simétrica** | $A = A^T$ |
| **Ortogonal** | $Q^{-1} = Q^T$ (rotación o reflexión) |
| **Similaridad** | $B = PAP^{-1}$; conserva valores propios |
| **Householder** | reflexión $I - \mathbf{u}\mathbf{u}^T/h$ para crear ceros |
| **QR** | $A = QR$, y se itera $RQ$ hasta triangular |
| **Tridiagonal** | solo 3 diagonales no nulas |
| **Hessenberg** | cero debajo de la primera subdiagonal |
| **Frobenius / compañera** | matriz cuyos valores propios son las raíces de un polinomio |
| **Diferencia hacia adelante** | $\Delta f_j = f_{j+1} - f_j$ |
| **Coeficientes de Markov** | coeficientes de $\binom{s}{k}$ en potencias de $s$ |
| **$\omega$** | frecuencia angular (rad/s); $\lambda = \omega^2$ |
| **Modo** | forma de vibrar asociada a un vector propio |
