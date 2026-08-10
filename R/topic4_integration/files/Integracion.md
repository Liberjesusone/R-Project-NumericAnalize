# Integración

**Carlos Echeverria**

**Universidad de "Los Andes".**

**Facultad de Ingeniería.**

**CESIMO (Centro de Simulación y Modelos)**

*March 7, 2022*

---

## Integración numérica

* Regla del trapecio

* Regla $1/3$ de Simpson

* Regla $3/8$ de Simpson

* Fórmula de Newton-Cotes

* Cuadratura de Gauss

---

## Regla del trapecio

*[Aquí va la gráfica: Gráfica de una función y = f(x) mostrando un trapecio sombreado bajo la curva entre los puntos a y b]*

Se obtiene al integrar la fórmula de interpolación lineal. Se escribe en la forma siguiente:

$$l = \int_{a}^{b} f(x) dx = \frac{b-a}{2}[f(a) + f(b)] + E \quad (1)$$

*[Aquí va la gráfica: Gráfica de y = f(x) con el área bajo la curva dividida en múltiples trapecios desde x_0 hasta x_N]*

La ecuación se puede extender a varios intervalos y se puede aplicar N veces, con una separación uniforme $h$

$$l = \int_{a}^{b} f(x) dx = \frac{b-a}{2}\left[f(a) + 2\sum_{j=1}^{N-1} f(a+jh) + f(b)\right] + E \quad (2)$$

donde $h = (b-a)/N$

### Ejemplo 4.1

*[Aquí va la gráfica: Figura E4.1 Un cuerpo de revolución]*

El cuerpo de revolución que se muestra en la figura E4.1 se obtiene al girar la curva dada por $y = 1 + (x/2)^2$, $0 \le x \le 2$, en torno al eje $x$. Calcule el volumen utilizando la regla extendida del trapecio con $N = 2, 4, 8, 16, 32, 64$ y $128$. El valor exacto es $I = 11.7286$. Evalúe el error para cada $N$.

**(Solución)**

El volumen está dado por

$$I = \int_{0}^{2} f(x) dx$$

donde

$$f(x) = \pi\left(1 + \left(\frac{x}{2}\right)^2\right)^2$$

A continuación aparecen los cálculos para $N = 2$ y $4$:

Para $N = 2$, $h = 2/2 = 1$

$$I = \frac{1}{2}[f(0) + 2f(1) + f(2)] = 0.5\pi[1 + 2(1.5625) + 4] = 12.7627$$

Para $N = 4$, $h = 2/4 = 0.5$

$$I = (0.5/2)[f(0) + 2f(0.5) + 2f(1) + 2f(1.5) + f(2)] = 11.9895$$

Las integraciones con los demás valores de $N$ se evalúan mediante el PROGRAMA 4-1. Los resultados se resumen en la tabla E4.2.

**Tabla E4.2**

| $N$ | $h$ | $I_h$ | $E_h$ |
| --- | --- | --- | --- |
| 2 | 1. | 12.7627 | 1.0341 |
| 4 | 0.5 | 11.9895 | -0.2609 |
| 8 | 0.25 | 11.7940 | -0.0654 |
| 16 | 0.125 | 11.7449 | -0.0163 |
| 32 | 0.0625 | 11.7326 | -0.0040 |
| 64 | 0.03125 | 11.7296 | -0.0010 |
| 128 | 0.015625 | 11.7288 | -0.0002 |

*Valor exacto: 11.7286*

El error de la regla del trapecio se define como

$$E = \int_{a}^{b} f(x) dx - \frac{b-a}{2}[f(a) + f(b)] \quad (3)$$

El desarrollo en serie de Taylor de $f(x)$, $f(a)$, $f(b)$ en torno a $\bar{x} = (a+b)/2$, con la hipótesis de que $f$ es analítica en $a \le x \le b$

$$f(x) = f(\bar{x}) + zf'(\bar{x}) + \frac{z^2}{2}f''(\bar{x}) + \dots$$

donde $z = x - \bar{x}$. El primer término de la ecuación

$$\int_{a}^{b} f(x) dx = \int_{-h/2}^{h/2} \left[f(\bar{x}) + zf'(\bar{x}) + \frac{z^2}{2}f''(\bar{x}) + \dots\right] dz \quad (4)$$

Al integrar obtenemos:

$$\int_{a}^{b} f(x) dx = hf(\bar{x}) + \frac{1}{24}h^3f''(\bar{x}) + \dots$$

Por otro lado el segundo término

$$\frac{b-a}{2}[f(a) + f(b)] = \frac{h}{2}\left[f(\bar{x}) - \frac{h}{2}f'(\bar{x}) + \frac{1}{2}\frac{h^2}{4}f''(\bar{x}) - \dots + f(\bar{x}) + \frac{h}{2}f'(\bar{x}) + \frac{1}{2}\frac{h^2}{4}f''(\bar{x}) + \dots\right]$$

y se obtiene

$$= hf(\bar{x}) + \frac{1}{8}h^3f''(\bar{x}) - \dots \quad (6)$$

$$E = \int_{a}^{b} f(x) dx - \frac{b-a}{2}[f(a) + f(b)] \simeq -\frac{1}{12}h^3f''(\bar{x}) \quad (7)$$

Supongamos que tenemos $N$ intervalos en $[a, b]$, entonces

$$E \simeq -\frac{1}{12}\frac{(b-a)^3}{N^3}\sum_{i=1}^{N} f''(\bar{x}_i) \quad (8)$$

ya que $h = (b-a)/N$ y $\bar{f}'' = N^{-1}\sum_{i=1}^{N} f''(\bar{x}_i)$. La ecuación queda

$$E \simeq -\frac{1}{12}(b-a)h^2\bar{f}''(\bar{x}_i) \quad (9)$$

válida si la función $f(x)$ es analítica en el intervalo.

Supongamos dos resultados $l_h$ y $l_{2h}$ y sus errores son

$$E_h \simeq Ch^2 \quad \text{y} \quad E_{2h} \simeq C(2h)^2 = C4h^2 \quad (10)$$

Por otra parte el valor exacto se puede escribir como

$$I = I_h + E_h = I_{2h} + E_{2h}$$

de lo cual obtenemos

$$E_h - E_{2h} = I_{2h} - I_h \quad (11)$$

Al despejar $C$

$$C = \frac{1}{3}h^{-2}(I_h - I_{2h}) \quad (12)$$

Así la primera ecuación es

$$E_h \simeq \frac{1}{3}(I_h - I_{2h})$$

De lo cual se obtiene una integral mas precisa:

$$I = I_h + E_h \simeq I_h + \frac{1}{3}(I_h - I_{2h}) \quad (13)$$

No es exacta pero el error es del orden $h^4$. Esta técnica se llama Integración de Romberg.

### Ejemplo 4.2

En el ejemplo 4.1, la regla extendida del trapecio da como resultado $I_{0.5} = 11.9895$ y $I_{0.25} = 11.7940$. Determine un valor más exacto utilizando la integración de Romberg.

**(Solución)**

Si definimos $h = 0.25$ en las ecuaciones (4.2.11) a la (4.2.14), el valor de $E_{0.25}$ dado por la ecuación (4.2.13) es

$$E_{0.25} \simeq \frac{1}{3}(11.7940 - 11.9895) = -0.0652$$

Por lo tanto, se obtiene el siguiente valor más exacto de $I$, dado por la ecuación (4.2.14)

$$I = I_{0.25} + E_{0.25} \simeq 11.7940 - 0.0652 = 11.7288$$

Este resultado coincide con el resultado para $N = 128$ (con $h = 0.0156$) en el ejemplo 4.1.

---

## Regla de $1/3$ de Simpson

*[Aquí va la gráfica: Gráfica mostrando el área bajo la curva aproximada con una parábola a través de tres puntos $f_0, f_1, f_2$]*

Se basa en la interpolación polinomial cuadrática. El polinomio de Newton hacia adelante ajustado a tres puntos, $x_0$, $x_1$, $x_2$.

$$l = \int_{a}^{b} f(x) dx = \frac{h^3}{3}[f(a) + 4f(\bar{x}) + f(b)] + E \quad (14)$$

donde $h = (b-a)/2$ y $\bar{x} = (a+b)/2$.

$$E \simeq -(h^5/90)f''''(\bar{x}) \quad (15)$$

La regla extendida es:

$$I = \int_{a}^{b} f(x) dx = \frac{h^3}{3}\left[f(a) + 4\sum_{i=1}^{N-1} f(a+ih) + 2\sum_{i=2}^{N-2} f(a+ih) + f(b)\right] + E \quad (16)$$

donde $h = (b-a)/N$. El error es

$$E \simeq -\frac{N}{2}\frac{h^5}{90}f''''(\bar{x}) = (a-b)\frac{h^4}{180}f''''(\bar{x}) \quad (17)$$

Para un dominio fijo $[a, b]$, el error es proporcional a $h^4$.

### Ejemplo 4.3

Repita el problema del ejemplo 4.1 utilizando la regla extendida de $1/3$ de Simpson con $N = 2, 4, 8, 16, 32$.

**(Solución)**

El tamaño del intervalo es $h = 2/N$. Los cálculos para $N = 2$ y $4$ son como sigue:

$N=2$

$$I = \frac{1}{3}[f(0) + 4f(1) + f(2)] = \frac{1}{3}\pi[1 + (4)(1.25^2) + 2^2] = 11.7809$$

$N=4$

$$I = \frac{0.5}{3}[f(0) + 4f(0.5) + 2f(1) + 4f(1.5) + f(2)] = \frac{0.5}{3}\pi[1 + 4(1.0625) + 2(1.25)^2 + 4(1.5625)^2 + 2^2] = 11.7318$$

Los cálculos para $N$ mayores se pueden realizar de manera análoga. Los resultados y evaluaciones del error se muestran a continuación.

| $N$ | $h$ | $I_h$ | $E_h$ |
| --- | --- | --- | --- |
| 2 | 1. | 11.7809 | -0.0523 |
| 4 | 0.5 | 11.7318 | -0.0032 |
| 8 | 0.25 | 11.7288 | -0.0002 |
| 16 | 0.125 | 11.7286 | 0 |
| 32 | 0.0625 | 11.7286 | 0 |
| 64 | 0.03125 | 11.7286 | 0 |

Al comparar los resultados anteriores con los del ejemplo 4.1 se puede ver que la regla extendida de Simpson es mucho más precisa que la regla extendida del trapecio, utilizando el mismo número de intervalos. Por ejemplo, la exactitud de la regla extendida del trapecio con $32$ intervalos es equivalente a la de la regla extendida de Simpson con tan sólo $4$ intervalos. El error de la regla extendida de Simpson es proporcional a $h^4$, por lo que es dos órdenes más grande que el de la regla extendida del trapecio.

Los desarrollos de Taylor para $f_0$ y $f_2$ en torno de $x_1$, o en forma equivalente, $\bar{x} = (a+b)/2$, se escriben como

$$f_0 = f_1 - hf_1' + \frac{1}{2}h^2f_1'' - \frac{1}{6}h^3f_1''' + \frac{1}{24}h^4f_1'''' - \dots$$

$$f_2 = f_1 + hf_1' + \frac{1}{2}h^2f_1'' + \frac{1}{6}h^3f_1''' + \frac{1}{24}h^4f_1'''' + \dots$$


Sustituyendolo en $l = \frac{1}{3}[f_0 + 4f_1 + f_2] + E$ obtenemos

$$I = 2hf_1 + \frac{1}{3}h^3f_1'' + \frac{1}{36}h^5f_1'''' + \dots + E \quad (18)$$

Por otro lado, el desarrollo de Taylor de $f(x)$ alrededor de $x_1$ es

$$f(x) = f_1 + zf_1' + \frac{1}{2}z^2f_1'' + \frac{1}{6}z^3f_1''' + \frac{1}{24}z^4f_1'''' - \dots$$

donde $x = x_1 + z$ o $z = x - x_1$.

La integración analítica de este desarrollo en $[a, b]$ da como resultado

$$\int_{a}^{b} f(x) dx = 2hf_1 + \frac{1}{3}h^3f_1'' + \frac{1}{60}h^5f_1'''' + \dots$$

Restamos y truncamos después del término principal, con lo que el error de la ecuación está dado aproximadamente por

$$E \simeq -\frac{1}{90}h^5f_1'''' \quad (19)$$

Una desventaja de la regla extendida de Simpson es que el número total de intervalos debe ser par. Por otro lado, la regla de $3/8$ de Simpson, se aplica únicamente a un número de intervalos que sea múltipo de tres. Por lo tanto, al combinar las reglas de $1/3$ y $3/8$, se puede considerar el caso tanto par como impar de intervalos.

---

## Regla de $3/8$ de Simpson

Para un dominio $[a, b]$ dividido en tres intervalos, se escribe como

$$l = \int_{a}^{b} f(x) dx = \frac{3}{8}h[f_0 + 3f_1 + 3f_2 + f_3] + E \quad (20)$$

donde $h = (b-a)/3$ y $f_i = f(a+ih)$. El error

$$E \simeq -\frac{3}{80}h^5f''''(\bar{x}) \quad (21)$$

con $\bar{x} = (a+b)/2$.

---

## Fórmulas de Newton-Cotes

Los métodos de integración numérica que se obtienen al integrar las formulas de interpolación de Newton reciben el nombre de formulas de Newton-Cotes. La regla del trapecio y las dos reglas de Simpson son casos de las fórmulas de Newton-Cotes, las cuales se dividen en fórmulas cerradas y abiertas.

Escribimos las fórmulas cerradas de Newton-Cotes en la forma:

$$\int_{a}^{b} f(x) dx = \alpha h[w_0f_0 + w_1f_1 + w_2f_2 + \dots + w_Nf_N] + E \quad (22)$$

donde $\alpha$ y $w$ son las constantes que aparecen en la tabla y $f_n = f(x_n)$, $x_n = a+nh$, $h = (b-a)/n$.

**Tabla 4.2 Constantes para las fórmulas cerradas de Newton-Cotes**

| $N$ | $\alpha$ | $w_i \quad (i = 0, 1, 2, \dots, N)$ | $E$ |
| --- | --- | --- | --- |
| $1$ | $1/2$ | $1, 1$ | $-\frac{1}{12}h^3f''$ |
| $2$ | $1/3$ | $1, 4, 1$ | $-\frac{1}{90}h^5f^{iv}$ |
| $3$ | $3/8$ | $1, 3, 3, 1$ | $-\frac{3}{80}h^5f^{iv}$ |
| $4$ | $2/45$ | $7, 32, 12, 32, 7$ | $-\frac{8}{945}h^7f^{vi}$ |
| $5$ | $5/288$ | $19, 75, 50, 50, 75, 19$ | $-\frac{275}{12096}h^7f^{vi}$ |
| $6$ | $1/140$ | $41, 216, 27, 272, 27, 216, 41$ | $-\frac{9}{1400}h^9f^{viii}$ |
| $7$ | $7/17280$ | $751, 3577, 1323, 2989, 2989, 1323, 3577, 751$ | $-\frac{8183}{518400}h^9f^{viii}$ |
| $8$ | $14/14175$ | $989, 5888, -928, 10946, -4540, 10946, -928, 5888, 989$ | $-\frac{2368}{467775}h^{11}f^{x}$ |
| $9$ | $9/89600$ | $2857, 15741, 1080, 19344, 5788, 5788, 19344, 1080, 15741, 2857$ | $-\frac{173}{14620}h^{11}f^{x}$ |
| $10$ | $5/299376$ | $16067, 106300, -48525, 272400, -260550, 272400, -48525, 106300, 16067$ | $-\frac{1346350}{326918592}h^{13}f^{xii}$ |

Recibe el nombre de fórmula cerrada, debido a que el dominio de integración está cerrado por el primer y último dato.

La integración se puede extender más allá de los puntos extremos (fórmulas abiertas). Dichas formulas se escriben como

$$\int_{a}^{b} f(x) dx = \alpha h[w_0f_0 + w_1f_1 + w_2f_2 + \dots + w_{N+2}f_{N+2}] + E \quad (23)$$


donde $h = (b-a)/(N+2)$. Las constantes $\alpha$ y $w$ se listan en la tabla, en donde $w_0$ y $w_{N+2}$ se igualan a cero debido a que corresponden a los extremos del dominio. Puesto que $w_0$ y $w_{N+2}$ se anulan, $f_0$ y $f_{N+2}$ son datos ficticios, que en realidad no son necesarios.

*[Aquí va la gráfica: Gráfica para las fórmulas abiertas de Newton-Cotes con puntos de interpolación interiores]*

Si comparamos una fórmula abierta con una cerrada utilizando, para el mismo número $N$ de datos, el error de la fórmula abierta es significativamente mayor que el de la fórmula cerrada. Por otro lado, se pueden utilizar las fórmulas abiertas cuando no se dispone de los valores de la función en los límites de integración.

**Tabla 4.3 Constantes para las fórmulas abiertas de Newton-Cotes**

| $N$ | $\alpha$ | $w_i \quad (i = 0, 1, \dots, N+2)$ | $E$ |
| --- | --- | --- | --- |
| $1$ | $3/2$ | $0, 1, 1, 0$ | $\frac{1}{4}h^3f''$ |
| $2$ | $4/3$ | $0, 2, -1, 2, 0$ | $\frac{28}{90}h^5f^{iv}$ |
| $3$ | $5/24$ | $0, 11, 1, 1, 11, 0$ | $\frac{95}{144}h^5f^{iv}$ |
| $4$ | $6/20$ | $0, 11, -14, 26, -14, 11, 0$ | $\frac{41}{140}h^7f^{vi}$ |
| $5$ | $7/1440$ | $0, 611, -453, 562, 562, -453, 611, 0$ | $\frac{5257}{8640}h^7f^{vi}$ |
| $6$ | $8/945$ | $0, 460, -954, 2196, -2459, 2196, -954, 460, 0$ | $\frac{3956}{14175}h^9f^{viii}$ |

### Ejemplo 4.4



*[Aquí va la gráfica: Figura E4.4 Aproximación lineal de una curva en coordenadas polares]*

La longitud de arco de una curva en coordenadas polares (véase la figura E4.4) está dada por

$$L = \int_{a}^{b} \sqrt{r^2 + \left(\frac{dr}{d\theta}\right)^2} d\theta$$


Calcular la longitud de arco de la curva dada por $r = 2(1+\cos\theta)$, $0 < \theta < \pi$, utilizando cada una de las fórmulas cerradas de integración de Newton-Cotes que aparecen en la tabla 4.2.

$$ds = \sqrt{(rd\theta)^2 + dr^2} = \sqrt{r^2 + (dr/d\theta)^2} d\theta$$

**(Solución)**


Hacemos $a = 0$ y $b = \pi$ y utilizamos el PROGRAMA 4-2. A continuación se muestran los resultados computacionales:

| Orden $N$ | Integral $L$ |
| --- | --- |
| 2 | 8.01823 |
| 3 | 8.00803 |
| 4 | 7.99993 |
| 5 | 7.99996 |
| 6 | 8.00000 |
| 7 | 8.00000 |
| 8 | 8.00000 |
| 9 | 8.00197 |
| 10 | 7.99201 |
| **Valor exacto** | **8.00000** |

Estos resultados ilustran el efecto de los errores por redondeo. Es decir, cuando $N$ crece, este resultado tiende el valor exacto, $8.00000$. Sin embargo, después de $N = 8$, el error vuelve a crecer de manera gradual. El crecimiento de los errores se atribuye a los errores de redondeo en las sumas y restas de los números demasiado grandes de las fórmulas.

---

## Cuadratura de Gauss



Utilizan puntos de Legendre (raíces de polinomios de Legendre). Las cuadraturas de Gauss no se pueden utilizar para integrar una función dada en forma de tabla con intervalos de separación uniforme debido a que los puntos de Legendre no están separados de esa manera; sin embargo, son más adecuadas para integrar funciones analíticas. La ventaja de las cuadraturas de Gauss es que su precisión es mayor que la de las formulas de Newton-Cotes.

El error de la regla del trapecio es proporcional a $f''$. Si se usa esta regla para integrar cada una de las funciones $f = 1, x, x^2, x^3, \dots$ entonces los resultados serán exactos para $f = 1$ y $f = x$, pero existirán errores para $x^2$ y las potencias superiores de $x$. El error de la regla de Simpson es proporcional a $f^{iv}$, por lo que es exacta si integramos $f = 1, x, x^2$ y $x^3$.

Sin embargo, el siguiente ejemplo muestra que una integración numérica con dos puntos se puede hacer exacta para el caso de los polinomios de orden tres, si se optimizan los valores $x$ de los datos.

### Ejemplo 4.5



La fórmula de integración con dos puntos se puede hacer exacta cuando se integra un polinomio de orden tres. Determine los puntos.

**(Solución)**


Consideremos


$$I = \int_{-1}^{1} f(x) dx \quad (A)$$


y escribamos una fórmula de integración con dos puntos, como sigue


$$I = w_1f(x_1) + w_2f(x_2) + E \quad (B)$$


donde $w_k$, $k=1, 2$, son pesos, $x_k$ son puntos indeterminados y $E$ es el término de error.

Ya que $w_k$ y $x_k$ son ambos indeterminados, requerimos que $E = 0$ (de modo que sea exacto) para $f(x) = 1, x, x^2$ y $x^3$. Introduciendo cada uno de los $f(x) = 1, x, x^2$ y $x^3$ en la ecuación (B) obtenemos cuatro ecuaciones:

$$2 = w_1 + w_2$$

$$0 = w_1x_1 + w_2x_2$$

$$\frac{2}{3} = w_1x_1^2 + w_2x_2^2$$

$$0 = w_1x_1^3 + w_2x_2^3$$


donde en el lado izquierdo aparecen los valores exactos.

Los límites de integración son $-1$ y $1$ y simétricos con respecto a $x = 0$, por lo que hacemos $x_2 = -x_1$ y requerimos que los puntos estén situados en forma simétrica. De la primera y segunda ecuación obtenemos


$$w_1 = w_2 = 1$$


Con estos valores, la cuarta ecuación se satisface automáticamente. La tercera ecuación es


$$\frac{1}{3} = x_1^2$$


de la cual se obtiene


$$x_1 = \frac{1}{\sqrt{3}} = 0.577350269$$


y


$$x_2 = -x_1 = -0.577350269$$

Con estos pesos y puntos, la ecuación (B) es exacta para un polinomio de orden menor o igual a tres. Aunque consideramos el intervalo $[-1, 1]$ por simplicidad, puede cambiarse por cualquier intervalo arbitrario mediante una transformación de coordenadas.

La fórmula general de la cuadratura de Gauss de orden $N$ es exacta cuando se integra un polinomio de orden menor o igual a $2N-1$.
Las cuadraturas de Gauss difieren en forma significativa de las formulas de Newton-Cotes ya que los $N$ puntos de la retícula (llamados puntos de Gauss) se obtienen mediante las raíces del polinomio de Legendre $P_N(x) = 0$ donde $P_N(x)$ es el polinomio de Legendre de orden $N$.

La cuadratura de Gauss que se extiende sobre el intervalo $[-1, 1]$ está dada por

$$\int_{-1}^{1} f(x) dx = \sum_{k=1}^{N} w_kf(x_k) \quad (24)$$


donde $N$ es el número de puntos de Gauss, los $w_i$ son los pesos y las $x_j$ son los puntos de Gauss dados en la tabla.

**Tabla 4.4 Puntos de Gauss y pesosª**

| $N$ | $\pm x_i$ | $w_i$ |
| --- | --- | --- |
| **$N = 2$** | 0.577350269 | 1.00000000 |
| **$N = 3$** | 0 | 0.888888889 |
|  | 0.774596669 | 0.555555556 |
| **$N = 4$** | 0.339981043 | 0.652145155 |
|  | 0.861136312 | 0.347854845 |
| **$N = 5$** | 0 | 0.568888889 |
|  | 0.538469310 | 0.478628670 |
|  | 0.906179846 | 0.236926885 |
| **$N = 6$** | 0.238619186 | 0.467913935 |
|  | 0.661209387 | 0.360761573 |
|  | 0.932469514 | 0.171324492 |
| **$N = 8$** | 0.183434642 | 0.362683783 |
|  | 0.525532410 | 0.313706646 |
|  | 0.796666478 | 0.222381034 |
|  | 0.960289857 | 0.101228536 |
| **$N = 10$** | 0.148874339 | 0.295524225 |
|  | 0.433395394 | 0.269266719 |
|  | 0.679409568 | 0.219086363 |
|  | 0.865063367 | 0.149451349 |
|  | 0.973906528 | 0.066671344 |

*ªVéase Abramowitz y Stegun para una tabla más completa.*

Por ejemplo, si $N = 4$, la ecuación es

$$\int_{-1}^{1} f(x) dx = 0.34785f(-0.86113) + 0.65214f(-0.33998) + 0.65214f(0.33998) + 0.34785f(0.86113) \quad (25)$$

La fórmula de integración de Gauss puede aplicarse a cualquier intervalo arbitrario $[a, b]$ con la transformación

$$x = \frac{2z-a-b}{b-a} \quad (26)$$


donde $z$ es la coordenada original en $a < z < b$ y $x$ es la coordenada normalizada en $-1 \le x \le 1$. La transformación de $x$ en $z$ es

$$z = \frac{(b-a)x + a + b}{2} \quad (27)$$

Por medio de esta transformación, la integral se puede escribir como

$$\int_{a}^{b} f(z) dz = \int_{-1}^{1} f(z)(dz/dx) dx = \frac{b-a}{2}\sum_{k=1}^{N} w_kf(z_k) \quad (28)$$


donde $dz/dx = (b-a)/2$. Los valores de $z_k$ se obtienen al sustituir $x_k$ en la ecuación (27) por los puntos de Gauss.

Por ejemplo, supongamos que $N = 2$, $a = 0$ y $b = 2$. Puesto que los puntos de Gauss $x_k$ para $N = 2$ en la coordenada normalizada $x$, $-1 \le x \le 1$, son $\pm0.57735$ (de la tabla), los puntos correspondientes en $z$ son

$$z_1 = \frac{1}{2}[(2-0)(-0.57735) + 0 + 2] = 0.42265$$

$$z_2 = \frac{1}{2}[(2-0)(0.57735) + 0 + 2] = 1.57735$$

La derivada es $dz/dx = (b-a)/2 = 1$. Por lo tanto, la cuadratura de Gauss se escribe como

$$\int_{a}^{b} f(z) dz = \int_{-1}^{1} f(z)(dz/dx) dx = [(1)f(0.42264) + (1)f(1.57735)] \quad (29)$$