# **Solución del Quiz: Interpolación Polinomial (Diferencias hacia adelante)**

## <font color="#00b0f0">Enunciado</font>
Genere 6 puntos igualmente espaciados para $f(x) = x\sin(x)$ en $\left[0,\,\dfrac{\pi}{2}\right]$ y construya la tabla de diferencias hacia adelante. Interpole $f\!\left(\dfrac{\pi}{4}\right)$ usando dos y cuatro nodos. Calcule el error respecto al valor real y analice las diferencias.

---
## <font color="#00b0f0">Generación de los 6 puntos</font>
6 puntos en $[0,\,\pi/2]$ implica 5 intervalos: $$h = \frac{\pi/2 - 0}{6-1} = \frac{\pi}{10} \approx 0.314159$$

| $i$ | $x_i$     | $x_i$ (decimal) |
| --- | --------- | --------------- |
| 0   | $0$       | $0.000000$      |
| 1   | $\pi/10$  | $0.314159$      |
| 2   | $2\pi/10$ | $0.628318$      |
| 3   | $3\pi/10$ | $0.942478$      |
| 4   | $4\pi/10$ | $1.256637$      |
| 5   | $\pi/2$   | $1.570796$      |

---
## <font color="#00b0f0">Evaluación de $f(x) = x\sin(x)$</font>
$$f_0 = 0 \cdot \sin(0) = 0.000000$$
$$f_1 = \tfrac{\pi}{10}\sin\!\left(\tfrac{\pi}{10}\right) = 0.314159 \times 0.309017 = 0.097081$$
$$f_2 = \tfrac{2\pi}{10}\sin\!\left(\tfrac{2\pi}{10}\right) = 0.628318 \times 0.587785 = 0.369316$$
$$f_3 = \tfrac{3\pi}{10}\sin\!\left(\tfrac{3\pi}{10}\right) = 0.942478 \times 0.809017 = 0.762482$$
$$f_4 = \tfrac{4\pi}{10}\sin\!\left(\tfrac{4\pi}{10}\right) = 1.256637 \times 0.951057 = 1.195135$$
$$f_5 = \tfrac{\pi}{2}\sin\!\left(\tfrac{\pi}{2}\right) = 1.570796 \times 1 = 1.570796$$

---
## <font color="#00b0f0">Tabla de diferencias hacia adelante</font>
La regla es siempre **siguiente menos actual**: $$\Delta f_i = f_{i+1} - f_i$$
**Orden 1:**
$$\Delta f_0 = 0.097081 - 0.000000 = 0.097081$$
$$\Delta f_1 = 0.369316 - 0.097081 = 0.272235$$
$$\Delta f_2 = 0.762482 - 0.369316 = 0.393166$$
$$\Delta f_3 = 1.195135 - 0.762482 = 0.432653$$
$$\Delta f_4 = 1.570796 - 1.195135 = 0.375661$$

**Orden 2** ($\Delta^2 f_i = \Delta f_{i+1} - \Delta f_i$):
$$\Delta^2 f_0 = 0.272235 - 0.097081 = 0.175154$$
$$\Delta^2 f_1 = 0.393166 - 0.272235 = 0.120931$$
$$\Delta^2 f_2 = 0.432653 - 0.393166 = 0.039487$$
$$\Delta^2 f_3 = 0.375661 - 0.432653 = -0.056992$$

**Orden 3:**
$$\Delta^3 f_0 = 0.120931 - 0.175154 = -0.054223$$
$$\Delta^3 f_1 = 0.039487 - 0.120931 = -0.081444$$
$$\Delta^3 f_2 = -0.056992 - 0.039487 = -0.096479$$

**Orden 4:**
$$\Delta^4 f_0 = -0.081444 - (-0.054223) = -0.027221$$
$$\Delta^4 f_1 = -0.096479 - (-0.081444) = -0.015035$$

**Orden 5:**
$$\Delta^5 f_0 = -0.015035 - (-0.027221) = 0.012186$$

### Tabla completa
En el examen solo realice las columnas necesarias para obtener el polinomio con 4 nodos

| $i$ |  $x_i$   |    $f_i$     |  $\Delta^1$  |  $\Delta^2$  |  $\Delta^3$   |  $\Delta^4$   |  $\Delta^5$  |
| :-: | :------: | :----------: | :----------: | :----------: | :-----------: | :-----------: | :----------: |
|  0  | 0.000000 | **0.000000** | **0.097081** | **0.175154** | **−0.054223** | **−0.027221** | **0.012186** |
|  1  | 0.314159 |   0.097081   |   0.272235   |   0.120931   |   −0.081444   |   −0.015035   |              |
|  2  | 0.628318 |   0.369316   |   0.393166   |   0.039487   |   −0.096479   |               |              |
|  3  | 0.942478 |   0.762482   |   0.432653   |  −0.056992   |               |               |              |
|  4  | 1.256637 |   1.195135   |   0.375661   |              |               |               |              |
|  5  | 1.570796 |   1.570796   |              |              |               |               |              |

Los coeficientes del polinomio son la **primera fila** (en negrita).

---
## <font color="#00b0f0">Cálculo de $s$</font>
Para interpolar en $x_A = \pi/4$, $s$ se calcula **una sola vez** respecto a $x_0$: $$s = \frac{x_A - x_0}{h} = \frac{\pi/4 - 0}{\pi/10} = \frac{\pi/4}{\pi/10} = \frac{10}{4} = 2.5$$
Este es uno de los puntos donde tuve más inflección ya qeu pensé que variaba $s$ según el nodo

---
## <font color="#00b0f0">Fórmula de Newton hacia adelante</font>
$$P_n(x_A) = f_0 + s\,\Delta f_0 + \frac{s(s-1)}{2!}\,\Delta^2 f_0 + \frac{s(s-1)(s-2)}{3!}\,\Delta^3 f_0 + \cdots$$

donde los factores numéricos con $s = 2.5$ son:

| Término          | Factor                          | Valor    |
| ---------------- | ------------------------------- | -------- |
| $s$              | $2.5$                           | $2.5$    |
| $s(s-1)/2!$      | $2.5 \times 1.5 / 2$            | $1.875$  |
| $s(s-1)(s-2)/3!$ | $2.5 \times 1.5 \times 0.5 / 6$ | $0.3125$ |

---
## <font color="#00b0f0">Interpolación con 2 nodos (grado 1)</font>
Usa $f_0$ y $\Delta f_0$:
$$P_1\!\left(\tfrac{\pi}{4}\right) = f_0 + s\,\Delta f_0 = 0 + 2.5 \times 0.097081 = \mathbf{0.242703}$$

---
## <font color="#00b0f0">Interpolación con 4 nodos (grado 3)</font>
Usa $f_0,\,\Delta f_0,\,\Delta^2 f_0,\,\Delta^3 f_0$:
$$P_3\!\left(\tfrac{\pi}{4}\right) = 0 + 2.5(0.097081) + 1.875(0.175154) + 0.3125(-0.054223)$$
$$= 0.242703 + 0.328414 - 0.016947 = \mathbf{0.554170}$$

---
## <font color="#00b0f0">Valor real y análisis del error</font>
$$f\!\left(\tfrac{\pi}{4}\right) = \tfrac{\pi}{4}\sin\!\left(\tfrac{\pi}{4}\right) = \tfrac{\pi}{4}\cdot\tfrac{\sqrt{2}}{2} = \tfrac{\pi\sqrt{2}}{8} \approx 0.555360$$

| Método          | $P({\pi}/{4})$ | Error absoluto | Error relativo |
| --------------- | -------------- | -------------- | -------------- |
| 2 nodos ($P_1$) | $0.242703$     | $0.312657$     | $56.3\%$       |
| 4 nodos ($P_3$) | $0.554170$     | $0.001190$     | $0.21\%$       |
| Exacto          | $0.555360$     | —              | —              |

### Análisis
Con 2 nodos el polinomio es lineal y usa únicamente $x_0 = 0$ y $x_1 = \pi/10$ para extrapolar hasta $\pi/4$, que se encuentra en $s = 2.5$, es decir, dos nodos y medio más allá del primer nodo. La extrapolación genera un error del $56\%$.

Con 4 nodos el polinomio cúbico utiliza los cuatro primeros nodos $\{0,\, \pi/10,\, 2\pi/10,\, 3\pi/10\}$. Ahora $\pi/4$ cae exactamente en el centro del intervalo cubierto ($s = 2.5$), dentro del rango de interpolación de los primeros tres nodos, y el polinomio captura la curvatura de $x\sin(x)$. El error cae a $0.21\%$ — una mejora de **más de 260 veces**.

Esto ilustra que:
1. La interpolación (dentro del intervalo de datos) es siempre más precisa que la extrapolación.
2. Aumentar el grado del polinomio incorpora órdenes superiores de diferencia que capturan la curvatura de la función.

---
## <font color="#00b0f0">Diferencias con mis respuestas</font>
- **1. Sentido de la resta en la tabla.**
	- Construí las diferencias en el sentido $f_i - f_{i+1}$ en vez de $f_{i+1} - f_i$. La convención estándar de Newton hacia adelante es siempre **siguiente menos actual** ($\Delta f_i = f_{i+1} - f_i$), que da valores positivos para una función creciente. El signo invertido afecta todos los coeficientes de la tabla y por ende el resultado de la interpolación.
- **2. La fórmula de evaluación usa un único $s$, no uno por nodo.**
	- Calculé $s_0, s_1, s_2$ distintos para cada término. La fórmula correcta usa **un solo valor** $s = (x_A - x_0)/h$ calculado respecto al primer nodo, y los factores sucesivos son $s-1,\, s-2,\ldots$ derivados de ese mismo $s$: $$P_n = f_0 + s\Delta f_0 + \frac{s(s-1)}{2!}\Delta^2 f_0 + \frac{s(s-1)(s-2)}{3!}\Delta^3 f_0 + \cdots$$
	- Confundí este esquema con el de diferencias divididas, donde los productos sí involucran $(x - x_i)$ explícitamente para cada nodo $i$.
- **3. Conclusión**
	- Debido a ambos errores y recordando el valor que presenté para el polinomio de 4 nodos, obtuvo un resultado similar a $4.6x10^{-3}$ es decir, bastante alejado de la aproximación que obtuvimos acá, al igual que en el polinomio de 2 nodos donde incluso obtuve un valor negativo. 
