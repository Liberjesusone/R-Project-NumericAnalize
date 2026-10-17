# Valores Propios

**Carlos Echeverria**  
*Universidad de "Los Andes"*  
*Facultad de Ingeniería*  
*CESIMO - Centro de Simulación y Modelos*  
*April 20, 2022*

---

## 1. Valores Propios de una Matriz

### Cálculo de valores propios de una matriz

Las ecuaciones lineales homogéneas se asocian con frecuencia a sistemas que presentan oscilaciones armónicas sin que haya fuerzas externas. En esta categoría están la vibración de una cuerda, una membrana u otros sistemas estructurales. Al estudiar la naturaleza o estabilidad dinámica de tales sistemas, es necesario resolver ecuaciones homogéneas.

Un conjunto en particular de ecuaciones homogéneas son:

$$3x - 2y + z = 0$$
$$x + y + 2z = 0$$
$$4x - y + 3z = 0$$

---

Si fijamos $x$ en un valor constante arbitrario $\beta$, tenemos:

$$-2y + z = -3\beta$$
$$y + 2z = -\beta$$
$$-y + 3z = -4\beta$$

Tiene solución si una de las ecuaciones es idéntica a otra, o bien si es combinación lineal de las demás.

---

Como el $\det(A) = 0$, el sistema de ecuaciones tiene solución. Donde $A$ es la matriz del sistema de ecuaciones y la solución es:

$$x = \beta, \quad y = \beta, \quad z = -\beta$$

Una forma estándar de un conjunto de ecuaciones homogéneas:

$$(a_{11} - \lambda)x_1 + a_{12}x_2 + a_{13}x_3 = 0$$
$$a_{21}x_1 + (a_{22} - \lambda)x_2 + a_{23}x_3 = 0$$
$$a_{31}x_1 + a_{32}x_2 + (a_{33} - \lambda)x_3 = 0$$

donde $\lambda$ es el valor característico o valor propio.

---

De forma equivalente:

$$Ax = \lambda x$$

donde

$$A = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}, \quad x = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}$$

Tiene una solución no trivial sólo si:

$$\det \begin{bmatrix} a_{11}-\lambda & a_{12} & a_{13} \\ a_{21} & a_{22}-\lambda & a_{23} \\ a_{31} & a_{32} & a_{33}-\lambda \end{bmatrix} = \det(A - \lambda I) = 0$$

---

Se puede definir $f(\lambda) \equiv \det(A - \lambda I)$, el cual es una función polinomial y su orden es el mismo que el de la matriz.

#### Ejemplo

<!-- ESPACIO PARA IMAGEN 1 -->
<!-- 
  [Pegar Imagen 1 Aquí]
  Descripción de la imagen: Esquema físico de un sistema mecánico compuesto por tres masas (m1, m2, m3) suspendidas verticalmente y unidas en serie por medio de resortes. Cada masa tiene indicado su desplazamiento vertical (y1, y2, y3).
-->

Sea un sistema vertical formado por masas y resortes. Al suponer que no hay fricción, las ecuaciones diferenciales para los desplazamientos de las masas son:

$$m_1 y_1''(t) = -(k_{01} + k_{12})y_1 + k_{12}y_2$$
$$m_2 y_2''(t) = k_{12}y_1 - (k_{12} + k_{23})y_2 + k_{23}y_3$$
$$m_3 y_3''(t) = k_{23}y_2 - (k_{23} + k_{34})y_3$$

Suponga $m_1 = m_2 = m_3 = m$.

---

La solución se puede escribir como:

$$y_j = \exp(i\omega t)f_j, \quad j = 1, 2, 3$$

donde $\omega$ es una velocidad angular no determinada, $f_j$ son las incógnitas e $i = \sqrt{-1}$.

$$-\omega^2 f_1 = -\left(\frac{1}{m}\right)(k_{01} + k_{12})f_1 + \left(\frac{1}{m}\right)k_{12}f_2$$
$$-\omega^2 f_2 = \left(\frac{1}{m}\right)k_{12}f_1 - \left(\frac{1}{m}\right)(k_{12} + k_{23})f_2 + \left(\frac{1}{m}\right)k_{23}f_3$$
$$-\omega^2 f_3 = \left(\frac{1}{m}\right)k_{23}f_2 - \left(\frac{1}{m}\right)(k_{23} + k_{34})f_3$$

---

donde $\lambda = \omega^2$

$$Af - \lambda f = 0$$

$$A = \begin{bmatrix} \frac{k_{01}+k_{12}}{m} & \frac{k_{12}}{m} & 0 \\ \frac{k_{21}}{m} & \frac{k_{12}+k_{23}}{m} & \frac{k_{23}}{m} \\ 0 & \frac{k_{32}}{m} & \frac{k_{23}+k_{34}}{m} \end{bmatrix}, \quad f = \begin{bmatrix} f_1 \\ f_2 \\ f_3 \end{bmatrix}$$

---

Debemos tomar en cuenta lo siguiente:

* Todos los valores propios de una matriz simétrica son reales. (Si todos los valores propios de una matriz simétrica son positivos, se dice que la matriz es definida positiva.)
* Una matriz no simétrica con entradas reales puede tener valores propios complejos, los cuales aparecen como parejas de complejos conjugados.

---

## 2. Transformación de Householder

Dada una matriz simétrica, ésta se puede transformar en una matriz tridiagonal.

* Transformación de una matriz a una tridiagonal:

$$A^{(1)} = A = \begin{bmatrix} x & x & x & \cdots \\ x & x & x & \cdots \\ x & x & x & \cdots \\ \vdots & \vdots & \vdots & \ddots \\ x & x & x & \cdots & x \end{bmatrix}$$

---

El primer paso para reducir $A^{(1)}$ a una matriz tridiagonal $A^{(2)}$:

$$A^{(2)} = \begin{bmatrix} x & x & 0 & \cdots & 0 \\ x & x & x & \cdots & x \\ 0 & x & x & \cdots & x \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & x & x & \cdots & x \end{bmatrix}$$

La transformación se lleva a cabo multiplicando por la izquierda y por la derecha por una matriz de transformación $P$:

$$A^{(2)} = PA^{(1)}P$$

---

$P$ se define con:

$$P = I - \frac{uu^T}{h}$$

$$u = \text{col}(0, a_{21} + G, a_{31}, a_{41}, \dots, a_{N1})$$

donde $a_{ij}$ son los elementos de $A^{(1)}$ y:

$$G = \left[ \sum_{i=2}^{N} (a_{i1})^2 \right]^{1/2} \text{signo}(a_{21})$$
$$h = G^2 + G a_{21}$$

---

donde $P$ tiene las siguientes características:

$$P = P^{-1}, \quad P^T = P, \quad PP = I$$

Repitiendo las transformaciones $N-2$ veces se consigue la matriz tridiagonal. En general, la transformación de $A^{(m)}$ en $A^{(m+1)}$ se escribe como:

$$A^{(m+1)} = P A^{(m)} P^{-1}$$

donde

$$P = I - \frac{uu^T}{h}$$

$$u = \text{col}(0, \dots, 0, a_{m+1, m} + G, a_{m+2, m}, \dots, a_{Nm})$$

donde $a_{ij}$ son los elementos de $A^{(m)}$ y:

---

$$G = \left[ \sum_{i=m+1}^{N} (a_{im})^2 \right]^{1/2} \text{signo}(a_{m+1, m})$$
$$h = G^2 + G a_{m+1, m}$$

La matriz tridiagonal se denota como:

$$M = \begin{bmatrix} a_1 & b_1 & \cdots & 0 \\ b_1 & a_2 & b_2 & \cdots \\ \vdots & b_2 & a_3 & b_3 \\ 0 & \cdots & b_N & a_N \end{bmatrix}$$

---

$$p_0(\lambda) = 1$$
$$p_1(\lambda) = a_1 - \lambda$$
$$p_2(\lambda) = (a_2 - \lambda)p_1(\lambda) - (b_1)^2$$
$$\vdots$$
$$p_i(\lambda) = (a_i - \lambda)p_{i-1}(\lambda) - (b_{i-1})^2 p_{i-2}(\lambda)$$

Para este caso se debe resolver $P_N(\lambda) = 0$.

* Raíces de $P_3(\lambda) = 0$: $\lambda_0, \lambda_1, \lambda_2$ (Valores propios de la tercera submatriz principal).
* Raíces de $P_4(\lambda) = 0$: $\lambda_0, \lambda_1, \lambda_2, \lambda_3$ (Valores propios de la cuarta submatriz principal).

---

#### Ejemplo:

$$A^{(1)} = A = \begin{bmatrix} 3 & 1 & 4 \\ 1 & 7 & 2 \\ 4 & 2 & 0 \end{bmatrix}$$

$$u_1 = 0$$
$$G = \sqrt{1^2 + 4^2} \text{signo}(1) = 4.1231$$
$$a_{21} = 1$$
$$h = (4.1231)^2 + (4.1231)(1) = 21.1231$$
$$u_2 = a_{21} + \sqrt{17} = 1 + 4.1231 = 5.1231$$
$$u_3 = a_{31} = 4$$

entonces $u = \text{col}[0, 5.1231, 4]$.

---

$$A^{(2)} = \left(I - \frac{uu^T}{h}\right) A^{(1)} \left(I - \frac{uu^T}{h}\right) = \begin{bmatrix} 3 & -4.1231 & 0 \\ -4.1231 & 1.3529 & 3.4118 \\ 0.0 & 3.4118 & 5.6471 \end{bmatrix}$$

El valor de $\lambda$ de $p_1(\lambda) = 0$ es $3.0$. Para $p_2(\lambda) = 0$ son $-2.0289$ y $6.3810$ y es equivalente a los autovalores de:

$$\begin{bmatrix} 3 & -4.1231 \\ -4.1231 & 1.3529 \end{bmatrix}$$

y para $p_3(\lambda) = 0$ tenemos $-2.8941$, $4.3861$ y $8.508$.

---

## 3. Método de las Potencias

Consideremos una matriz $A$. Los valores y vectores propios satisfacen:

$$A u_i = \lambda_i u_i$$

donde $\lambda_i$ es el i-ésimo valor propio y $u_i$ es el i-ésimo vector propio. Supondremos que el máximo valor propio es real y aislado y que los valores propios están numerados en orden creciente:

$$|\lambda_1| \le |\lambda_2| \le \cdots \le |\lambda_{N-1}| < |\lambda_N|$$

El método comienza con una estimación inicial del vector propio, $u^{(0)}$, que puede ser cualquier vector no nulo. La primera aproximación iterativa es:

$$u^{(1)} = A u^{(0)}$$

---

y las siguientes:

$$u^{(k+1)} = \frac{1}{\lambda^{(k)}} A u^{(k)}$$

con

$$\lambda^{(k)} = \frac{(u^{(k)}, u^{(k)})}{(u^{(k)}, \frac{1}{\lambda^{(k-1)}} u^{(k-1)})}$$

donde $(a, b)$ denota el producto escalar y $k$ el número de la iteración.

---

#### Ejemplo

$$A = \begin{bmatrix} 3 & 1 & 4 \\ 1 & 7 & 2 \\ 4 & 2 & 0 \end{bmatrix}, \quad u^{(0)} = \text{col}(1, 1, 1)$$

**Tabla 7.3 Solución iterativa:**

| $k$ | $\lambda^{(k)}$ | Elementos de $u^{(k)}$ |
| :---: | :---: | :---: |
| 1 | 8.333333 | 0.9600000, 0.7200000, 1.200000 |
| 2 | 8.453462 | 0.8233314, 0.7381592, 1.277583 |
| 3 | 8.492605 | 0.7889469, 1.323826, 0.6886570 |
| 4 | 8.503810 | 0.7579303, 0.6824518, 1.344461 |
| 5 | 8.506876 | 0.7462269, 0.6724730, 1.355852 |
| 10 | 8.507998 | 0.7319254, 1.366952, 0.6655880 |
| 15 | 8.508000 | 0.7314215, 0.6653084, 1.367357 |
| 18 | 8.508000 | 0.7314049, 1.367370, 0.6652992 |
| 19 | 8.508000 | 0.7314036, 0.6652986, 1.367371 |
| 20 | 8.508000 | 0.7314029, 0.6652982, 1.367372 |

---

### Método de las potencias inversas
Para el menor valor propio:

$$|\lambda_1| < |\lambda_2| \le \cdots \le |\lambda_{N-1}| \le |\lambda_N|$$

El primer paso:

$$A u^{(1)} = u^{(0)}$$

y las iteraciones:

$$u^{(k+1)} = \frac{1}{\lambda^{(k)}} A u^{(k)}$$

con

$$\lambda^{(k)} = \frac{(u^{(k)}, \lambda^{(k-1)} u^{(k-1)})}{(u^{(k)}, u^{(k)})}$$

---

## 4. Iteración QR

Si se denota la matriz inicial $A = A_0 = Q_0 R_0$, donde $Q_0$ es una matriz ortogonal y $R_0$ es una matriz triangular superior. La transformación de similaridad se escribe como:

$$A_1 = Q_0^{-1} A_0 Q_0$$

o, en forma equivalente:

$$A_1 = R_0 Q_0$$

Para mejorar la eficiencia de los cálculos, se hacen dos modificaciones. La primera es:

$$A_k - \alpha_k I = Q_k R_k$$
$$A_{k+1} = R_k Q_k + \alpha_k I \quad (1)$$

---

La segunda se refiere a aplicar el método de Householder al caso de una matriz no simétrica:

$$\begin{bmatrix} x & x & x & x & x \\ x & x & x & x & x \\ 0 & x & x & x & x \\ 0 & 0 & x & x & x \\ 0 & 0 & 0 & x & x \\ 0 & 0 & 0 & \vdots & x \\ 0 & 0 & 0 & 0 & x \end{bmatrix}$$

que se llama forma superior de Hessenberg.

---

La iteración QR transforma cualquier matriz real a una matriz triangular superior de bloques, como:

$$\begin{bmatrix} D & x & x & x & x \\ 0 & D & x & x & x \\ 0 & 0 & D & x & x \\ 0 & 0 & 0 & D & x \\ 0 & 0 & 0 & 0 & D \end{bmatrix}$$

donde $x$ y $D$ representan los elementos no nulos. La matriz completa tiene forma de matriz triangular superior por bloques.

---

Los valores propios de la anterior matriz son iguales a los de una matriz diagonal por bloques y se obtiene haciendo cero todos los $x$:

$$\begin{bmatrix} D & 0 & 0 & 0 & 0 \\ 0 & D & 0 & 0 & 0 \\ 0 & 0 & D & 0 & 0 \\ 0 & 0 & 0 & D & 0 \\ 0 & 0 & 0 & 0 & D \end{bmatrix}$$

Los valores propios son iguales a los de los bloques diagonales.

---

Una aplicación importante de la iteración QR es el cálculo de raíces de un polinomio. La matriz de Frobenius se define como:

$$P = \begin{bmatrix} -p_{N-1} & -p_{N-2} & -p_{N-3} & \cdots & -p_0 \\ 1 & 0 & 0 & \cdots & 0 \\ 0 & 1 & 0 & \cdots & 0 \\ 0 & \vdots & 1 & \cdots & 0 \\ 0 & 0 & 0 & 1 & 0 \end{bmatrix}$$

su función característica es:

$$f(\lambda) = \det \begin{bmatrix} -p_{N-1}-\lambda & -p_{N-2} & -p_{N-3} & \cdots & -p_0 \\ 1 & -\lambda & 0 & \cdots & 0 \\ 0 & 1 & -\lambda & \cdots & 0 \\ 0 & \vdots & 1 & \ddots & 0 \\ 0 & 0 & 0 & 1 & -\lambda \end{bmatrix}$$

---

Desarrollando en series de potencias:

$$f(\lambda) = |P - \lambda I| = (-1)^N [\lambda^N + p_{N-1}\lambda^{N-1} + \cdots + p_0]$$

Por otro lado, cualquier ecuación polinomial:

$$a_N x^N + \cdots + a_2 x^2 + a_1 x + a_0 = 0$$

se puede expresar en la forma:

$$p_0 = a_0 / a_N$$
$$p_1 = a_1 / a_N$$
$$p_2 = a_2 / a_N$$
$$\vdots$$
$$p_{N-1} = a_{N-1} / a_N$$

Por lo tanto se puede utilizar la iteración QR para hallar sus autovalores.

---

#### Ejemplo:
Calcule las raíces de la siguiente ecuación:

$$f(x) = x^5 - 0.2x^4 + 7x^3 + x^2 - 3.5x + 2 = 0$$

La matriz de Frobenius:

$$M = \begin{bmatrix} 0.2 & -7 & -1 & 3.5 & -2 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \end{bmatrix}$$

los autovalores son: $-0.9085$, $0.1403 \pm 2.7314i$ y $0.4139 \pm 0.3506i$.