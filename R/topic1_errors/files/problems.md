# PROBLEMAS

**1.1)** Si se usan 8 bits para representar los enteros positivos y negativos en complemento a dos, ¿cuál es el entero positivo más grande y el negativo más pequeño (en magnitud) en decimal?

**1.2)** Se tienen dos números binarios de 16 bits en complemento a dos:

**a)** Binario: `1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1`
&nbsp;&nbsp;&nbsp;&nbsp;(bit no: `0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15`).

**b)** Binario: `1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0`
&nbsp;&nbsp;&nbsp;&nbsp;(bit no: `0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15`)

Determine los valores decimales de los dos números binarios.

**1.3)** Hallar el épsilon de la máquina para una IBM PC en WATFOR-77.

**1.4)** Repita el problema anterior con la computadora mainframe a la que tenga acceso.

**1.5)** Evalúe
$$ \exp(x) - 1 $$

para $x = 0.0001$, aplicando el desarrollo de Taylor para $\exp(x)$. Use los primeros tres términos.

**1.6)** Desarrolle las siguientes funciones en serie de Maclaurin.

$$ 1/(1 + x^2) $$

$$ \tan(x) $$

$$ 1/(1 - x) $$

$$ \ln(1 + x) $$

**1.7)** Muestre que el desarrollo de Taylor de $\ln[(1 + x)/(1 - x)]$ alrededor de $x = 0$ es

$$ 2 \sum_{n=1}^{\infty} \frac{x^{2n-1}}{2n-1} $$

**1.8)** Por medio del desarrollo de Maclaurin de $e^x$ y $e^{-x}$, obtenga el desarrollo de Maclaurin de $\text{senh} (x)$ y $\text{cosh} (x)$, donde

$$ \text{senh} (x) = \frac{1}{2}(e^x - e^{-x}) $$

$$ \text{cosh} (x) = \frac{1}{2}(e^x + e^{-x}) $$

**1.9 a)** Si la siguiente función se escribe en un programa, ¿en cuál rango de $x$ aparecerá un desborde o una división entre cero originados por el error de redondeo?

$$ f(x) = \frac{1}{1 - \tanh(x)} $$

Suponga que el número positivo más pequeño es $3 \times 10^{-39}$ y el épsilon de la máquina es $1.2 \times 10^{-7}$.

**b)** Reescriba la ecuación de tal forma que no se necesite restar.
