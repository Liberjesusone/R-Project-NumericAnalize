# Álgebra Lineal

**Autor:** Carlos Echeverria  
**Institución:** Universidad de "Los Andes", Facultad de Ingeniería  
**Fecha:** 1 de abril de 2022  

---

## 1. Álgebra Lineal Numérica

El álgebra lineal es fundamental, tanto para el análisis científico como para los métodos numéricos, y no podríamos hacer mucho sin tener un conocimiento básico de ella.

---

## 2. Eliminación de Gauss y Gauss-Jordan

Dado un conjunto de $N$ ecuaciones:

$$\begin{aligned}
a_{1,1}x_{1} + a_{1,2}x_{2} + a_{1,3}x_{3} + \dots + a_{1,N}x_{N} &= y_{1} \\
a_{2,1}x_{1} + a_{2,2}x_{2} + a_{2,3}x_{3} + \dots + a_{2,N}x_{N} &= y_{2} \\
&\vdots \\
a_{N,1}x_{1} + a_{N,2}x_{2} + a_{N,3}x_{3} + \dots + a_{N,N}x_{N} &= y_{N}
\end{aligned} \tag{1}$$

donde $a_{i,j}$ son los coeficientes, $x_{j}$ son las incógnitas y los $y_{i}$ son los términos independientes. Si al menos uno de $\{y_{j}\}$ es diferente de cero, el conjunto es no homogéneo. Considerando que el conjunto de ecuaciones tiene solución única, la eliminación de Gauss consiste en:

1. **La eliminación hacia adelante**
2. **La sustitución hacia atrás**

### Eliminación hacia adelante

Comenzamos multiplicando la primera ecuación por $\frac{a_{2,1}}{a_{1,1}}$, luego se le resta a la segunda ecuación y se sustituye en la misma, y así sucesivamente con el resto de las ecuaciones:

$$\begin{aligned}
a_{1,1}x_{1} + a_{1,2}x_{2} + a_{1,3}x_{3} + \dots + a_{1,N}x_{N} &= y_{1} \\
a'_{2,2}x_{2} + a'_{2,3}x_{3} + \dots + a'_{2,N}x_{N} &= y'_{2} \\
&\vdots \\
a'_{N,2}x_{2} + a'_{N,3}x_{3} + \dots + a'_{N,N}x_{N} &= y'_{N}
\end{aligned} \tag{2}$$

donde:

$$a'_{i,j} = a_{i,j} - \left(\frac{a_{i,1}}{a_{1,1}}\right)a_{1,j}$$

Repitiendo el procedimiento en todas las restantes ecuaciones:

$$\begin{aligned}
a_{1,1}x_{1} + a_{1,2}x_{2} + a_{1,3}x_{3} + \dots + a_{1,N}x_{N} &= y_{1} \\
a'_{2,2}x_{2} + a'_{2,3}x_{3} + \dots + a'_{2,N}x_{N} &= y'_{2} \\
a''_{3,3}x_{3} + \dots + a''_{3,N}x_{N} &= y''_{3} \\
&\vdots \\
a_{N,N}^{(N-1)}x_{N} &= y_{N}^{(N-1)}
\end{aligned} \tag{3}$$

### Sustitución hacia atrás

Los términos principales reciben el nombre de **pivotes**. El procedimiento de sustitución hacia atrás comienza con la última ecuación:

$$\begin{aligned}
x_{N} &= \frac{y_{N}^{(N-1)}}{a_{N,N}^{(N-1)}} \\
x_{N-1} &= \frac{y_{N-1}^{(N-2)} - a_{N-1,N}^{(N-2)}x_{N}}{a_{N-1,N-1}^{(N-2)}} \\
&\vdots \\
x_{1} &= \frac{y_{1} - \sum_{j=2}^{N} a_{1,j}x_{j}}{a_{1,1}}
\end{aligned} \tag{4}$$

### Representación Matricial (Arreglo Aumentado)

Para mayor facilidad, es preferible escribir las ecuaciones en un arreglo:

$$\begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} & \dots & a_{1,N-1} & a_{1,N} & y_{1} \\
a_{2,1} & a_{2,2} & a_{2,3} & \dots & a_{2,N-1} & a_{2,N} & y_{2} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
a_{N,1} & a_{N,2} & a_{N,3} & \dots & a_{N,N-1} & a_{N,N} & y_{N}
\end{bmatrix}$$

Y con las etapas de la eliminación:

$$\begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} & \dots & a_{1,N-1} & a_{1,N} & y_{1} \\
0 & a'_{2,2} & a'_{2,3} & \dots & a'_{2,N-1} & a'_{2,N} & y'_{2} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \dots & a_{N-1,N-1}^{(N-2)} & a_{N-1,N}^{(N-1)} & y_{N-1}^{(N-2)} \\
0 & 0 & 0 & \dots & 0 & a_{N,N}^{(N-1)} & y_{N}^{(N-1)}
\end{bmatrix} \tag{5}$$

---

### Ejemplo de Eliminación de Gauss

Dado el sistema:

$$\begin{aligned}
2x_{1} + x_{2} - 3x_{3} &= -1 \\
-x_{1} + 3x_{2} + 2x_{3} &= 12 \\
3x_{1} + x_{2} - 3x_{3} &= 0
\end{aligned}$$

Expresión en arreglo:

$$\left[ \begin{array}{ccc|c}
2 & 1 & -3 & -1 \\
-1 & 3 & 2 & 12 \\
3 & 1 & -3 & 0
\end{array} \right]$$

**Procedimiento:**
* $\text{Fila } 2 - \text{Fila } 1 \cdot (a_{2,1}/a_{1,1})$
* $\text{Fila } 3 - \text{Fila } 1 \cdot (a_{3,1}/a_{1,1})$

$$\left[ \begin{array}{ccc|c}
2 & 1 & -3 & -1 \\
0 & 7/2 & 1/2 & 23/2 \\
0 & -1/2 & 3/2 & 3/2
\end{array} \right] \longrightarrow \left[ \begin{array}{ccc|c}
2 & 1 & -3 & -1 \\
0 & 7/2 & 1/2 & 23/2 \\
0 & 0 & 11/7 & 22/7
\end{array} \right]$$

Sustituyendo hacia atrás obtenemos:
$$x_{3} = 2, \quad x_{2} = 3, \quad x_{1} = 1$$

---

### Método de Gauss-Jordan

La eliminación de Gauss-Jordan hacia atrás difiere del anterior, ya que el proceso es en sentido inverso. Partiendo de la ecuación en la forma reducida (5), convierte en $1$ los coeficientes en la posición de pivoteo y elimina a los demás:

$$\begin{bmatrix} 0 & 0 & 0 & \dots & 0 & 1 & \bar{y}_{N} \end{bmatrix}$$

donde $\bar{y}_{N} = y_{N}^{(N-1)} / a_{N,N}^{(N-1)}$.

Los $N$-ésimos coeficientes de cada renglón, excepto el último, se eliminan restando el último renglón:

$$\begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} & \dots & a_{1,N-1} & 0 & \bar{y}_{1} \\
0 & a'_{2,2} & a'_{2,3} & \dots & a'_{2,N-1} & 0 & \bar{y}_{2} \\
0 & 0 & a''_{3,3} & \dots & a'_{3,N-1} & 0 & \bar{y}_{3} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \dots & a_{N-1,N-1}^{(N-2)} & 0 & \bar{y}_{N-1} \\
0 & 0 & 0 & \dots & 0 & 1 & \bar{y}_{N}^{(N-1)}
\end{bmatrix}$$

donde $\bar{y}_{i} = y_{i}^{(i-1)} - a_{i,N}^{(i-1)}\bar{y}_{N}$. Siguiendo un procedimiento análogo:

$$\begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} & \dots & 0 & 0 & \bar{y}_{1} \\
0 & a'_{2,2} & a'_{2,3} & \dots & 0 & 0 & \bar{y}_{2} \\
0 & 0 & a''_{3,3} & \dots & 0 & 0 & \bar{y}_{3} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \dots & 1 & 0 & \bar{y}'_{N-1} \\
0 & 0 & 0 & \dots & 0 & 1 & \bar{y}_{N}
\end{bmatrix}$$

Finalmente:

$$\begin{bmatrix}
1 & 0 & 0 & \dots & 0 & 0 & \bar{y}_{1}^{(N-1)} \\
0 & 1 & 0 & \dots & 0 & 0 & \bar{y}_{2}^{(N-2)} \\
0 & 0 & 1 & \dots & 0 & 0 & \bar{y}_{3}^{(N-3)} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \dots & 1 & 0 & \bar{y}'_{N-1} \\
0 & 0 & 0 & \dots & 0 & 1 & \bar{y}_{N}
\end{bmatrix}$$

Es de notar que todos los coeficientes son iguales a cero, excepto los pivotes. Cada renglón se interpreta como:

$$x_{i} = \bar{y}_{i}^{(N-i)} \tag{6}$$

---

### Ejemplo de Gauss-Jordan

Resuelva el ejemplo anterior con el método Gauss-Jordan:

$$\begin{aligned}
\text{Fila 1:} \quad &\begin{bmatrix} 2 & 1 & -3 & -1 \end{bmatrix} \\
\text{Fila 2:} \quad &\begin{bmatrix} 0 & 7/2 & 1/2 & 23/2 \end{bmatrix} \\
\text{Fila 3:} \quad &\begin{bmatrix} 0 & 0 & 11/7 & 22/7 \end{bmatrix}
\end{aligned}$$

* $\text{Fila 3} \div a_{3,3} \implies \begin{bmatrix} 0 & 0 & 1 & 2 \end{bmatrix} \equiv \text{Fila C}$
* Repitiendo:
  * $\text{Fila 1} - a_{1,3}\text{Fila C} = \begin{bmatrix} 2 & 1 & 0 & 5 \end{bmatrix}$
  * $\text{Fila 2} - a_{2,3}\text{Fila C} = \begin{bmatrix} 0 & 7 & 0 & 21 \end{bmatrix}$
  * $\text{Fila 3} \div a_{3,3} = \begin{bmatrix} 0 & 0 & 1 & 2 \end{bmatrix}$

Resultando en:

$$\begin{bmatrix}
2 & 0 & 0 & 2 \\
0 & 1 & 0 & 3 \\
0 & 0 & 1 & 2
\end{bmatrix} \longrightarrow \begin{bmatrix}
1 & 0 & 0 & 1 \\
0 & 1 & 0 & 3 \\
0 & 0 & 1 & 2
\end{bmatrix}$$

---

## 3. Pivoteo y Eliminación de Gauss

El método de Gauss no funciona si el primer coeficiente de la primera fila es igual a cero o si un coeficiente de la diagonal se anula en el proceso de la solución. El **pivoteo** se usa para solventar estos inconvenientes.

Consideremos el arreglo:

$$\left[ \begin{array}{ccc|c}
0 & 10 & 1 & 2 \\
1 & 3 & -1 & 6 \\
2 & 4 & 1 & 5
\end{array} \right] \xrightarrow{\text{Reordenando}} \left[ \begin{array}{ccc|c}
2 & 4 & 1 & 5 \\
1 & 3 & -1 & 6 \\
0 & 10 & 1 & 2
\end{array} \right]$$

Por eliminación:

$$\left[ \begin{array}{ccc|c}
2 & 4 & 1 & 5 \\
0 & 1 & -3/2 & 7/2 \\
0 & 10 & 1 & 2
\end{array} \right] \longrightarrow \left[ \begin{array}{ccc|c}
2 & 4 & 1 & 5 \\
0 & 10 & 1 & 2 \\
0 & 1 & -3/2 & 7/2
\end{array} \right]$$

Por eliminación posterior:

$$\left[ \begin{array}{ccc|c}
2 & 4 & 1 & 5 \\
0 & 10 & 1 & 2 \\
0 & 0 & -16/10 & 33/10
\end{array} \right]$$

La sustitución hacia atrás da:
* $x_{3} = -2.0625$
* $x_{2} = \frac{2 - x_{3}}{10} = 0.4062$
* $x_{1} = \frac{5 - 4x_{2} - x_{3}}{2} = 2.7187$

Las eliminaciones de Gauss-Jordan llevan a:

$$\begin{bmatrix}
1 & 0 & 0 & 2.7187 \\
0 & 1 & 0 & 0.4062 \\
0 & 0 & 1 & -2.0625
\end{bmatrix}$$

---

## 4. Problemas sin Solución Única

Ejemplos:

a)
$$\begin{aligned}
-x + y &= 1 \\
-2x + 2y &= 2
\end{aligned}$$

b)
$$\begin{aligned}
-x + y &= 1 \\
-x + y &= 0
\end{aligned}$$

c)
$$\begin{aligned}
-x + y &= 1 \\
x + 2y &= -2 \\
2x - y &= 0
\end{aligned}$$

---

### Representación Gráfica

<!-- IMAGEN 1 -->
![Interpretación Geométrica de Sistemas de Ecuaciones](espacio_para_imagen_1)

> **Descripción de la Imagen 1:**
> Tres gráficas en el plano cartesiano $xy$ representando los tipos de soluciones geométricas para sistemas de ecuaciones lineales de $2 \times 2$:
> 1. **Número infinito de soluciones:** Dos ecuaciones representan exactamente la misma recta que se superpone (coincidentes).
> 2. **Sin solución:** Dos rectas paralelas con la misma pendiente y distinta intersección, por lo que jamás se cruzan.
> 3. **Sin solución (sistema sobredeterminado):** Tres rectas en el plano que se intersectan en puntos distintos de a pares, pero no poseen un único punto común de intersección de las tres a la vez.

---

## 5. Vectores y Matrices

### Definiciones

Una Matriz $A$:

$$A = \begin{bmatrix}
2 & 1 & 0 \\
-1 & 3 & 1 \\
1 & 1 & 1
\end{bmatrix} = \begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix} = A[a_{ij}]$$

Una matriz rectangular de dimensión $m \times n$ (ejemplo $2 \times 3$):

$$A = \begin{bmatrix}
5 & 9 & 2 \\
3 & 0 & 4
\end{bmatrix}$$

### Operaciones Básicas

* **Suma:** $A + B = C \implies c_{ij} = a_{ij} + b_{ij}$
* **Resta:** $A - B = C \implies c_{ij} = a_{ij} - b_{ij}$
* **Multiplicación:** $AB = C \implies c_{ij} = \sum_{k=1}^{N} a_{ik}b_{kj}$

### Propiedades de las Operaciones

* El producto $AB \neq BA$ en general. Si $AB = BA$, se dice que las matrices **conmutan**.
* Si $A$ y $B$ son rectangulares, el producto solo existe si $A$ es de dimensión $m \times n$ y $B$ es de dimensión $n \times k$.
* **División:** $B^{-1}A = C \implies A = BC$. La división es más restrictiva, puesto que $B^{-1}$ solo puede existir si $B$ es cuadrada.

### Vectores

Un vector columna:

$$\mathbf{x} = \begin{bmatrix} x_{1} \\ x_{2} \\ x_{3} \end{bmatrix}, \quad \mathbf{y} = \begin{bmatrix} y_{1} \\ y_{2} \\ y_{3} \end{bmatrix}$$

Un vector renglón:

$$\mathbf{a} = \begin{bmatrix} a_{1} & a_{2} & a_{3} & a_{4} \end{bmatrix}$$

* **Suma de vectores:** $\mathbf{x} + \mathbf{y} = \mathbf{z} \implies x_{i} + y_{i} = z_{i}$
* **Resta de vectores:** $\mathbf{x} - \mathbf{y} = \mathbf{z} \implies x_{i} - y_{i} = z_{i}$
* **Multiplicación Matriz-Vector:** $A\mathbf{x} = \mathbf{y} \implies y_{i} = \sum_{k=1}^{N} a_{ik} x_{k}$

---

### Ejemplo de Operaciones Matriciales

Sean:

$$A = \begin{bmatrix}
1 & 2 & 4 \\
3 & 1 & 2 \\
4 & 1 & 3
\end{bmatrix}, \quad
B = \begin{bmatrix}
7 & 3 & 1 \\
2 & 3 & 5 \\
8 & 1 & 6
\end{bmatrix}, \quad
\mathbf{x} = \begin{bmatrix} 1 \\ 4 \\ 2 \end{bmatrix}, \quad
\mathbf{y} = \begin{bmatrix} 3 \\ 9 \\ 4 \end{bmatrix}$$

**Calcular $A+B$, $A-B$, $AB$ y $BA$:**

$$A+B = \begin{bmatrix}
1+7 & 2+3 & 4+1 \\
3+2 & 1+3 & 2+5 \\
4+8 & 1+1 & 3+6
\end{bmatrix} = \begin{bmatrix}
8 & 5 & 5 \\
5 & 4 & 7 \\
12 & 2 & 9
\end{bmatrix}$$

$$A-B = \begin{bmatrix}
1-7 & 2-3 & 4-1 \\
3-2 & 1-3 & 2-5 \\
4-8 & 1-1 & 3-6
\end{bmatrix} = \begin{bmatrix}
-6 & -1 & 3 \\
1 & -2 & -3 \\
-4 & 0 & -3
\end{bmatrix}$$

$$AB = \begin{bmatrix}
1 & 2 & 4 \\
3 & 1 & 2 \\
4 & 1 & 3
\end{bmatrix}
\begin{bmatrix}
7 & 3 & 1 \\
2 & 3 & 5 \\
8 & 1 & 6
\end{bmatrix} = \begin{bmatrix}
1\cdot7+2\cdot2+4\cdot8 & 1\cdot3+2\cdot3+4\cdot1 & 1\cdot1+2\cdot5+4\cdot6 \\
3\cdot7+1\cdot2+2\cdot8 & 3\cdot3+1\cdot3+2\cdot1 & 3\cdot1+1\cdot5+2\cdot6 \\
4\cdot7+1\cdot2+3\cdot8 & 4\cdot3+1\cdot3+3\cdot1 & 4\cdot1+1\cdot5+3\cdot6
\end{bmatrix} = \begin{bmatrix}
43 & 13 & 35 \\
39 & 14 & 20 \\
54 & 18 & 27
\end{bmatrix}$$

$$BA = \begin{bmatrix}
7 & 3 & 1 \\
2 & 3 & 5 \\
8 & 1 & 6
\end{bmatrix}
\begin{bmatrix}
1 & 2 & 4 \\
3 & 1 & 2 \\
4 & 1 & 3
\end{bmatrix} = \begin{bmatrix}
20 & 18 & 37 \\
31 & 12 & 29 \\
35 & 23 & 52
\end{bmatrix}$$

---

### Tipos de Matrices Especiales

* **Matriz Nula:**

$$A = \begin{bmatrix}
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{bmatrix}$$

* **Matriz Identidad:**

$$I = \begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

* **Matriz Transpuesta:**

$$A = \begin{bmatrix} 2 & 3 \\ 0 & 5 \end{bmatrix} \implies A^t = \begin{bmatrix} 2 & 0 \\ 3 & 5 \end{bmatrix}$$

$$B = \begin{bmatrix} 5 & 0 \\ 2 & 7 \\ 1 & 2 \end{bmatrix} \implies B^t = \begin{bmatrix} 5 & 2 & 1 \\ 0 & 7 & 2 \end{bmatrix}$$

* **Vector Transpuesto:**

$$\mathbf{x} = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} \implies \mathbf{x}^t = \begin{bmatrix} x_1 & x_2 & x_3 \end{bmatrix}$$

* **Matriz Inversa:** Satisface $A A^{-1} = A^{-1} A = I$.

* **Vector Nulo:** $\mathbf{x} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}$

* **Vectores Unitarios:**

$$\mathbf{u} = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}, \quad \mathbf{v} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}, \quad \mathbf{w} = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$$

---

## 6. Inversión de una Matriz

Consideremos la ecuación lineal:

$$A\mathbf{x} = \mathbf{y} \tag{7}$$

donde $A$ es una matriz cuadrada. Supongamos que el pivoteo no es necesario:

$$G A \mathbf{x} = G \mathbf{y} \tag{8}$$

Si escogemos $G$ como la inversa de $A$ ($G = A^{-1}$):

$$\mathbf{x} = A^{-1}\mathbf{y} \tag{9}$$

La eliminación de Gauss-Jordan es equivalente a la multiplicación a la izquierda por $G = A^{-1}$. Escribimos $A$ e $I$ en forma de arreglo aumentado $[A \mid I]$:

$$\left[ \begin{array}{ccc|ccc}
a_{1,1} & a_{1,2} & a_{1,3} & 1 & 0 & 0 \\
a_{2,1} & a_{2,2} & a_{2,3} & 0 & 1 & 0 \\
a_{3,1} & a_{3,2} & a_{3,3} & 0 & 0 & 1
\end{array} \right] \longrightarrow \left[ \begin{array}{ccc|ccc}
1 & 0 & 0 & g_{1,1} & g_{1,2} & g_{1,3} \\
0 & 1 & 0 & g_{2,1} & g_{2,2} & g_{2,3} \\
0 & 0 & 1 & g_{3,1} & g_{3,2} & g_{3,3}
\end{array} \right]$$

---

### Ejemplo de Inversión de Matriz

Calcule la inversa de la matriz $A$:

$$A = \begin{bmatrix}
2 & 1 & -3 \\
-1 & 3 & 2 \\
3 & 1 & -3
\end{bmatrix}$$

**Paso 1:** Planteamiento en arreglo aumentado $[A \mid I]$:

$$\left[ \begin{array}{ccc|ccc}
2 & 1 & -3 & 1 & 0 & 0 \\
-1 & 3 & 2 & 0 & 1 & 0 \\
3 & 1 & -3 & 0 & 0 & 1
\end{array} \right]$$

**Paso 2:** Eliminación de los elementos debajo de la primera columna:

$$\left[ \begin{array}{ccc|ccc}
2 & 1 & -3 & 1 & 0 & 0 \\
0 & 3.5 & 0.5 & 0.5 & 1 & 0 \\
0 & -0.5 & 1.5 & -1.5 & 0 & 1
\end{array} \right]$$

**Paso 3:** Eliminación del elemento debajo del segundo pivote:

$$\left[ \begin{array}{ccc|ccc}
2 & 1 & -3 & 1 & 0 & 0 \\
0 & 3.5 & 0.5 & 0.5 & 1 & 0 \\
0 & 0 & 1.5714 & -1.4285 & 0.14285 & 1
\end{array} \right]$$

Dividiendo la tercera fila por $1.5714$:

$$\left[ \begin{array}{ccc|ccc}
2 & 1 & -3 & 1 & 0 & 0 \\
0 & 3.5 & 0.5 & 0.5 & 1 & 0 \\
0 & 0 & 1 & -0.90909 & 0.90909 & 0.63636
\end{array} \right]$$

**Paso 4:** Eliminación hacia arriba (Gauss-Jordan). Dividir el segundo renglón entre $3.5$ y restar el último renglón multiplicado por $0.5/3.5$:

$$\left[ \begin{array}{ccc|ccc}
2 & 1 & -3 & 1 & 0 & 0 \\
0 & 1 & 0 & 0.27272 & 0.27272 & -0.09090 \\
0 & 0 & 1 & -0.90909 & 0.90909 & 0.63636
\end{array} \right]$$

Seguidamente, normalizando y eliminando la fila superior:

$$\left[ \begin{array}{ccc|ccc}
1 & 0 & 0 & -1 & 0 & 1 \\
0 & 1 & 0 & 0.27272 & 0.27272 & -0.09090 \\
0 & 0 & 1 & -0.90909 & 0.90909 & 0.63636
\end{array} \right]$$

### Verificación del resultado:

$$A A^{-1} = \begin{bmatrix}
2 & 1 & -3 \\
-1 & 3 & 2 \\
3 & 1 & -3
\end{bmatrix}
\begin{bmatrix}
-1 & 0 & 1 \\
0.27272 & 0.27272 & -0.09090 \\
-0.90909 & 0.90909 & 0.63636
\end{bmatrix} = \begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

$$A^{-1} A = \begin{bmatrix}
-1 & 0 & 1 \\
0.27272 & 0.27272 & -0.09090 \\
-0.90909 & 0.90909 & 0.63636
\end{bmatrix}
\begin{bmatrix}
2 & 1 & -3 \\
-1 & 3 & 2 \\
3 & 1 & -3
\end{bmatrix} = \begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

---

## 7. Descomposición LU

Es una transformación de una matriz $A$ como producto de dos matrices:

$$A = LU \tag{10}$$

donde $L$ es una matriz triangular inferior con unos en la diagonal y $U$ es una matriz triangular superior.

$$\begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} \\
a_{2,1} & a_{2,2} & a_{2,3} \\
a_{3,1} & a_{3,2} & a_{3,3}
\end{bmatrix} = \begin{bmatrix}
1 & 0 & 0 \\
l_{2,1} & 1 & 0 \\
l_{3,1} & l_{3,2} & 1
\end{bmatrix}
\begin{bmatrix}
u_{1,1} & u_{1,2} & u_{1,3} \\
0 & u_{2,2} & u_{2,3} \\
0 & 0 & u_{3,3}
\end{bmatrix}$$

### Esquema General de la Descomposición LU

a) **Primer renglón de $U$** ($j = 1, \dots, N$):  
$$u_{1,j} = a_{1,j}, \quad j = 1, \dots, N$$

b) **Primera columna de $L$** ($i = 2, \dots, N$):  
$$l_{i,1} = \frac{a_{i,1}}{u_{1,1}}, \quad i = 2, \dots, N$$

c) **Segundo renglón de $U$** ($j = 2, \dots, N$):  
$$u_{2,j} = a_{2,j} - l_{2,1}u_{1,j}, \quad j = 2, \dots, N$$

d) **Segunda columna de $L$** ($i = 3, \dots, N$):  
$$l_{i,2} = \frac{a_{i,2} - l_{i,1}u_{1,2}}{u_{2,2}}, \quad i = 3, \dots, N$$

e) **El $n$-ésimo renglón de $U$** ($j = n, \dots, N$):  
$$u_{n,j} = a_{n,j} - \sum_{k=1}^{n-1} l_{n,k}u_{k,j}, \quad j = n, \dots, N$$

f) **La $n$-ésima columna de $L$** ($i = n+1, \dots, N$):  
$$l_{i,n} = \frac{a_{i,n} - \sum_{k=1}^{n-1} l_{i,k}u_{k,n}}{u_{n,n}}, \quad i = n+1, \dots, N$$

---

### Ejemplo de Descomposición LU

Dada la matriz $A$:

$$A = \begin{bmatrix}
2 & 1 & -3 \\
-1 & 3 & 2 \\
3 & 1 & -3
\end{bmatrix}$$

**Cálculos:**
* $u_{1,1} = 2, \quad u_{1,2} = 1, \quad u_{1,3} = -3$
* $l_{2,1} = -0.5, \quad l_{3,1} = 1.5$
* $u_{2,2} = 3 - (-0.5)(1) = 3.5$
* $u_{2,3} = 2 - (-0.5)(-3) = 0.5$
* $l_{3,2} = \frac{1 - (1.5)(1)}{3.5} = -0.142857$
* $u_{3,3} = -3 - (1.5)(-3) - (-0.142857)(0.5) = 1.57142$

**Resultado:**

$$L = \begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
1.5 & -0.1429 & 1
\end{bmatrix}, \quad
U = \begin{bmatrix}
2 & 1 & -3 \\
0 & 3.5 & 0.5 \\
0 & 0 & 1.5714
\end{bmatrix}$$

Matriz compacta combinada:

$$\begin{bmatrix}
2 & 1 & -3 \\
-0.5 & 3.5 & 0.5 \\
1.5 & -0.1429 & 1.5714
\end{bmatrix} = \begin{bmatrix}
u_{1,1} & u_{1,2} & u_{1,3} \\
l_{2,1} & u_{2,2} & u_{2,3} \\
l_{3,1} & l_{3,2} & u_{3,3}
\end{bmatrix}$$

---

### Solución de Sistemas con Descomposición LU

Dada la ecuación:

$$A\mathbf{x} = LU\mathbf{x} = \mathbf{y} \tag{11}$$

Sea $U\mathbf{x} = \mathbf{z}$, lo que implica $L\mathbf{z} = \mathbf{y}$: $\tag{12}$

$$\begin{bmatrix}
1 & 0 & 0 \\
l_{2,1} & 1 & 0 \\
l_{3,1} & l_{3,2} & 1
\end{bmatrix}
\begin{bmatrix} z_1 \\ z_2 \\ z_3 \end{bmatrix} =
\begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} \tag{13}$$

Sustitución hacia adelante para hallar $\mathbf{z}$:

$$\begin{aligned}
z_{1} &= y_{1} \\
z_{2} &= y_{2} - z_{1}l_{2,1} \\
z_{3} &= y_{3} - z_{1}l_{3,1} - z_{2}l_{3,2}
\end{aligned}$$

Luego, se resuelve $U\mathbf{x} = \mathbf{z}$:

$$\begin{bmatrix}
u_{1,1} & u_{1,2} & u_{1,3} \\
0 & u_{2,2} & u_{2,3} \\
0 & 0 & u_{3,3}
\end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} =
\begin{bmatrix} z_1 \\ z_2 \\ z_3 \end{bmatrix} \tag{14}$$

Sustitución hacia atrás para hallar $\mathbf{x}$:

$$\begin{aligned}
x_{3} &= \frac{z_{3}}{u_{3,3}} \\
x_{2} &= \frac{z_{2} - u_{2,3}x_{3}}{u_{2,2}} \\
x_{1} &= \frac{z_{1} - u_{1,2}x_{2} - u_{1,3}x_{3}}{u_{1,1}}
\end{aligned}$$

---

### Algoritmo General para LU

$$\begin{aligned}
z_{1} &= y_{1} \\
z_{i} &= y_{i} - \sum_{j=1}^{i-1} l_{i,j}z_{j}, \quad i = 2, 3, \dots, N \\
x_{N} &= \frac{z_{N}}{u_{N,N}} \\
x_{i} &= \frac{z_{i} - \sum_{j=i+1}^{N} u_{i,j}x_{j}}{u_{i,i}}, \quad i = N-1, N-2, \dots, 1
\end{aligned}$$

---

## 8. Determinantes

El determinante es un número característico asociado con toda matriz cuadrada:

$$\det(A) = \sum (\pm) a_{i1}a_{j2}a_{k3}\dots a_{rN} \tag{15}$$

donde la suma se realiza sobre todas las permutaciones de los primeros subíndices de $a$, asignando el signo $(+)$ si la permutación es par y $(-)$ si la permutación es impar.

### Determinante $2 \times 2$

$$\det(A) = \det \begin{bmatrix} a_{1,1} & a_{1,2} \\ a_{2,1} & a_{2,2} \end{bmatrix} = a_{1,1}a_{2,2} - a_{1,2}a_{2,1}$$

---

### Determinante $3 \times 3$ (Regla de Sarrus)

$$\begin{aligned}
\det(A) &= \det \begin{bmatrix}
a_{1,1} & a_{1,2} & a_{1,3} \\
a_{2,1} & a_{2,2} & a_{2,3} \\
a_{3,1} & a_{3,2} & a_{3,3}
\end{bmatrix} \\
&= a_{1,1}a_{2,2}a_{3,3} + a_{2,1}a_{3,2}a_{1,3} + a_{3,1}a_{1,2}a_{2,3} \\
&\quad - a_{1,1}a_{3,2}a_{2,3} - a_{2,1}a_{1,2}a_{3,3} - a_{3,1}a_{2,2}a_{1,3}
\end{aligned}$$

Una matriz de orden $N$ posee $N!$ permutaciones.

<!-- IMAGEN 2 -->
![Regla de Sarrus Esquema](espacio_para_imagen_2)

> **Descripción de la Imagen 2:**
> Esquema gráfico con el patrón de cálculo de la Regla de Sarrus para una matriz de $3 \times 3$. Muestra las entradas $a_{11}, a_{12}, a_{13}$, $a_{21}, a_{22}, a_{23}$, y $a_{31}, a_{32}, a_{33}$ unidas por líneas discontinuas diagonales:
> * Las tres diagonales descendentes hacia la derecha representan los productos con signo positivo $(+)$.
> * Las tres diagonales ascendentes hacia la derecha (o descendentes hacia la izquierda) representan los productos con signo negativo $(-)$.

---

### Cálculo de Determinantes mediante Descomposición LU

Otra forma eficiente de calcular el determinante es utilizar la eliminación de Gauss o la descomposición $LU$:

* **Regla 1:** $\det(BC) = \det(B)\det(C)$
* **Regla 2:** El $\det(M)$ es igual al producto de todos los elementos de la diagonal principal de $M$, si $M$ es una matriz triangular (superior o inferior).

Usando la descomposición $LU$:

$$\det(A) = \det(LU) = \det(L)\det(U) = \det(U) \tag{16}$$

puesto que $\det(L) = 1$. Por lo tanto:

$$\det(U) = \prod_{k=1}^{N} u_{k,k} \tag{17}$$

---

### Ejemplo

Dada la matriz $A$:

$$A = \begin{bmatrix}
2 & 1 & -3 \\
-1 & 3 & 2 \\
3 & 1 & -3
\end{bmatrix}$$

Sus matrices $L$ y $U$ son:

$$L = \begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
1.5 & -0.1429 & 1
\end{bmatrix}, \quad
U = \begin{bmatrix}
2 & 1 & -3 \\
0 & 3.5 & 0.5 \\
0 & 0 & 1.5714
\end{bmatrix}$$

$$\det(L) = 1$$

$$\det(U) = (2)(3.5)(1.5714) = 11$$

Por lo tanto, $\det(A) = 11$.

---

## 9. Problemas Mal Condicionados

Ciertos problemas tienen solución, aunque sus soluciones se vuelven muy imprecisas debido a severos errores de redondeo. Tienen las siguientes características:

* Un ligero cambio en los coeficientes provoca cambios significativos en la solución.
* Los elementos de la diagonal de la matriz tienden a ser menores que los elementos que no pertenecen a la diagonal.
* El producto $\det(A) \cdot \det(A^{-1})$ difiere en forma significativa de $1$.
* El resultado de $(A^{-1})^{-1}$ es muy distinto de $A$.
* $AA^{-1}$ difiere en alto grado de la matriz identidad $I$.
* $A^{-1}(A^{-1})^{-1}$ difiere más de la matriz identidad que lo que difiere $AA^{-1}$.

---

### Ejemplo: Matriz de Hilbert $4 \times 4$

Sea la matriz de Hilbert de dimensiones $4 \times 4$:

$$A = \begin{bmatrix}
1 & 1/2 & 1/3 & 1/4 \\
1/2 & 1/3 & 1/4 & 1/5 \\
1/3 & 1/4 & 1/5 & 1/6 \\
1/4 & 1/5 & 1/6 & 1/7
\end{bmatrix}$$

**Ejercicio propuesto:**
Calcular:
a) $A^{-1}(A^{-1})^{-1}$  
b) $\det(A) \cdot \det(A^{-1})$