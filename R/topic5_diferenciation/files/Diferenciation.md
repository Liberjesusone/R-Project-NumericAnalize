# Diferenciación 
Carlos Echeverria 
Universidad de "Los Andes". 
Facultad de Ingeniería. 
March 23, 2022 

## Contenido 
* Definición de diferencial 
* Uso del desarrollo de Taylor 
* Algoritmo genérico para una aproximación por diferencias 
* Uso de los operadores de diferencia 
* Aproximación de derivadas parciales por diferencias 

## Diferenciación numérica 
También llamada aproximación por diferencias, se utiliza para evaluar las derivadas de una función por medio de sus valores dados en los puntos de una retícula. 

**Diferencias hacia adelante** 
$$f'(x_0) \approx \frac{f(x_0 + h) - f(x_0)}{h}$$ 

**Diferencias hacia atrás** 
$$f'(x_0) \approx \frac{f(x_0) - f(x_0 - h)}{h}$$ 

**Diferencias centrales** 
$$f'(x_0) \approx \frac{f(x_0 + h) - f(x_0 - h)}{2h}$$ 

<!-- [MARCA: Insertar imagen de las gráficas de Hacia adelante, Hacia atrás y Central aquí] --> 

## Uso del desarrollo de Taylor 
Cuando una función se representa numéricamente en puntos discretos, esta se aproxima mediante la interpolación.  Para la obtención de la aproximación por diferencias se utiliza el desarrollo de Taylor que es equivalente a la diferenciación por interpolación.  Otra particularidad es que para una derivada de orden $p$, el mínimo número de datos necesario para obtener una aproximación por diferencias es $p + 1$. 

Comencemos con definir $f'_i = f'(x_i)$, $f_i = f(x_i)$ y $f_{i+1} = f(x_{i+1})$.  El desarrollo de Taylor: 
$$f_{i+1} = f_i + hf'_i + \frac{h^2}{2}f''_i + \frac{h^3}{6}f'''_i + \frac{h^4}{24}f^{(4)}_i + \dots$$ 

Despejando: 
$$f'_i = \frac{f_{i+1} - f_i}{h} - \frac{h}{2}f''_i - \frac{h^2}{6}f'''_i - \dots$$ 

El error por truncamiento: 
$$f'_i = \frac{f_{i+1} - f_i}{h} + O(h)$$ 
donde $O(h) = -\frac{h}{2}f''_i$.  

La aproximación por diferencias hacia atrás de la primera derivada: 
$$f_{i-1} = f_i - hf'_i + \frac{h^2}{2}f''_i - \frac{h^3}{6}f'''_i + \frac{h^4}{24}f^{(4)}_i - \dots$$ 

Al despejar: 
$$f'_i = \frac{f_i - f_{i-1}}{h} + O(h)$$ 
donde $O(h) = \frac{h}{2}f''_i$.  

Si restamos: 
$$f_{i+1} - f_{i-1} = 2hf'_i + \frac{h^3}{3}f'''_i + \dots$$ 
$$f'_i = \frac{f_{i+1} - f_{i-1}}{2h} - \frac{h^2}{6}f'''_i + \dots$$ 

La aproximación por diferencias centrales se expresa como: 
$$f'_i = \frac{f_{i+1} - f_{i-1}}{2h} + O(h^2)$$ 
donde $O(h^2) = -\frac{h^2}{6}f'''_i$. 

Usemos: 
$$f_{i+1} = f_i + hf'_i + \frac{h^2}{2}f''_i + \frac{h^3}{6}f'''_i + \frac{h^4}{24}f^{(4)}_i + \dots$$ 
$$f_{i+2} = f_i + 2hf'_i + 4\frac{h^2}{2}f''_i + 8\frac{h^3}{6}f'''_i + 16\frac{h^4}{24}f^{(4)}_i + \dots$$ 
para eliminar los términos de la segunda derivada: 
$$4f_{i+1} - f_{i+2} = 3f_i + 2hf'_i - \frac{2}{3}h^3f'''_i + \dots$$ 

Al despejar $f'_i$: 
$$f'_i = \frac{-f_{i+2} + 4f_{i+1} - 3f_i}{2h} + O(h^2)$$ 
donde $O(h^2) = -\frac{1}{3}h^2f'''_i$. 

La ecuación se llama aproximación por diferencias hacia adelante con tres puntos para $f'_i$ y el error es del mismo orden que el de la aproximación por diferencias centrales. 

Análogamente, la aproximación por diferencias hacia atrás con tres puntos: 
$$f'_i = \frac{3f_i - 4f_{i-1} + f_{i-2}}{2h} + O(h^2)$$ 
donde $O(h^2) = \frac{1}{3}h^2f'''_i$. 

<!-- [MARCA: Insertar imagen de Ejemplo 5.1 aquí] --> 

Conviene observar que los errores de las dos primeras aproximaciones decrecen en proporción con $h$, mientras que los errores de las últimas tres aproximaciones decrecen en proporción con $h^2$.  Es claro que a razón de reducción del error se vuelve más rápida cuando el orden de precisión es mayor. 

<!-- [MARCA: Insertar imagen de la tabla de resultados de Ejemplo 5.1 aquí] --> 

Para obtener una segunda derivada: 
$$f_{i+1} + f_{i-1} = 2f_i + h^2f''_i + \frac{1}{12}h^4f^{(4)}_i + \dots$$ 

Restando $2f_i$: 
$$f_{i+1} - 2f_i + f_{i-1} = h^2f''_i + \frac{1}{12}h^4f^{(4)}_i + \dots$$ 

Entonces: 
$$f''_i = \frac{f_{i+1} - 2f_i + f_{i-1}}{h^2} + O(h^2)$$   
donde $O(h^2) = -\frac{1}{12}h^2f^{(4)}_i$.  y se llama aproximación por diferencias centrales de $f''_i$. 

También se puede usar $f_i$, $f_{i-1}$ y $f_{i-2}$: 
$$f_{i-2} - 2f_{i-1} = -f_i + h^2f''_i - h^3f'''_i + \dots$$ 

Al despejar: 
$$f''_i = \frac{f_{i-2} - 2f_{i-1} + f_i}{h^2} + O(h)$$   
donde $O(h) = hf'''_i$.  
se llama aproximación por diferencias hacia atrás de $f''_i$. 

<!-- [MARCA: Insertar imágenes de la Tabla 5.2 (Aproximaciones por diferencias) partes a-i aquí] --> 

## Algoritmo para una aproximación por diferencias 
El objetivo de esta sección es describir un algoritmo genérico para obtener una aproximación por diferencias de una derivada de orden dado, utilizando un conjunto específico de puntos en una retícula. 

Supongamos que el número total de puntos en la retícula es $L$ y que los puntos de la retícula son $i = \alpha, \beta, \dots, \lambda$.  Supongamos que $L > p + 1$, donde $p$ es el orden de la derivada por aproximar.  Las abscisas de los puntos de la retícula son $x_i = \alpha h, \beta h, \dots, \lambda h$. 

La aproximación por diferencias de la $p$-ésima derivada de $f(x)$, utilizando estos puntos de la retícula, se puede escribir en la forma: 
$$f^{(p)}_0 = \frac{a_\alpha f_\alpha + a_\beta f_\beta + \dots + a_\lambda f_\lambda}{h^p} + E$$   
donde $a_\alpha$ hasta $a_\lambda$ son $L$ coeficientes indeterminados y 
$$E = c_1 h^{L-p} f^{(L)} + c_2 h^{L-p+1} f^{(L+1)}$$   

<!-- [MARCA: Insertar imagen de la Figura 5.2 (Ilustración de los puntos de la retícula) aquí] --> 

La esencia del algoritmo es sustituir los desarrollos de Taylor de $f$, y calcular los coeficientes indeterminados de forma que el término del error se minimice o, en forma equivalente, que el orden de $E$ sea el máximo orden posible. 

Para simplificar la explicación posterior, supongamos que $p = 1$, $L = 3$, $\alpha = 0$, $\beta = 1$ y $\lambda = 2$.  Entonces: 
$$f'_0 = \frac{a_0 f_0 + a_1 f_1 + a_2 f_2}{h} + E$$   
donde $x_0 = 0$, $x_1 = h$ y $x_2 = 2h$.  Sustituimos los desarrollos de Taylor de $f_1$ y $f_2$ alrededor de $x = 0$: 
$$f'_0 = \frac{1}{h} \left[ a_0 f_0 + a_1 \left( f_0 + hf'_0 + \frac{h^2}{2}f''_0 + \frac{h^3}{6}f'''_0 + \dots \right) + a_2 \left( f_0 + 2hf'_0 + \frac{4h^2}{2}f''_0 + \frac{8h^3}{6}f'''_0 + \dots \right) \right] + E$$ 

Reagrupamos términos: 
$$f'_0 = \frac{f_0[a_0 + a_1 + a_2]}{h} + f'_0[0 + a_1 + 2a_2] + \frac{f''_0 h[0 + a_1 + 4a_2]}{2} + \frac{f'''_0 h^2[0 + a_1 + 8a_2]}{6} + \frac{f^{(4)}_0 h^3[0 + a_1 + 16a_2]}{24} + \dots + E$$ 

Hay tres coeficientes indeterminados, los cuales se pueden definir mediante tres condiciones.  Para minimizar el error, hacemos los coeficientes de $f_0$, $f'_0$ y $f''_0$ iguales a 0, 1 y 0 respectivamente: 
$$a_0 + a_1 + a_2 = 0$$ 
$$0 + a_1 + 2a_2 = 1$$   
$$0 + a_1 + 4a_2 = 0$$ 

Al resolver las ecuaciones anteriores, vemos que los valores de los tres coeficientes indeterminados son: $a_0 = -3/2$, $a_1 = 2$ y $a_2 = -1/2$. 

Los términos no nulos de orden superior constituyen el error; 
$$E = -f'''_0 \frac{h^2 [0 + a_1 + 8a_2]}{6} - f^{(4)}_0 \frac{h^3 [0 + a_1 + 16a_2]}{24} + \dots$$   

Al comparar con la ecuación  : 
$$c_1 = -\frac{1}{6}(a_1 + 8a_2) = \frac{1}{3}$$ 
$$c_2 = -\frac{1}{24}(a_1 + 16a_2) = \frac{1}{4}$$ 
$$E = \frac{1}{3}h^2f'''_0$$ 

El resultado final es: 
$$f'_0 = \frac{1}{h} \left( -\frac{3}{2}f_0 + 2f_1 - \frac{1}{2}f_2 \right) + E = \frac{-3f_0 + 4f_1 - f_2}{2h} + E$$   

En términos más generales, con $L$ datos, podemos definir $L$ coeficientes indeterminados, fijando correctamente los primeros $L$ términos del desarrollo de Taylor.  Así, el término del error es proporcional al $(L+1)$-ésimo término o, en forma equivalente, a la $L$-ésima derivada si su coeficiente no se anula.  Si es igual a cero, el término del error es del orden inmediato superior. 

Este algoritmo funciona incluso cuando los índices $\alpha, \beta, \dots$ no son enteros.  Esto quiere decir que la aproximación por diferencias en una retícula con separación no uniforme se puede obtener mediante el mismo algoritmo. 

## Uso de los operadores de diferencias 
Definimos a continuación tres operadores de diferencias: 

Operador de diferencias hacia adelante: $\Delta$ 
$$\Delta f_i = f_{i+1} - f_i$$   

Operador de diferencias hacia atrás: $\nabla$ 
$$\nabla f_i = f_i - f_{i-1}$$   

Operador diferencial central: $\delta$ 
$$\delta f_i = f_{i+\frac{1}{2}} - f_{i-\frac{1}{2}}$$  o bien 
$$\delta f_{i+\frac{1}{2}} = f_{i+1} - f_i$$   
donde $f_{i+\frac{1}{2}} = f(x_i + h/2)$. 

Los operadores de diferencias de orden superior se pueden escribir como potencias de los operadores de diferencias anteriores: por ejemplo, $\Delta^n$, $\nabla^n$ y $\delta^n$ son operadores de diferencias de orden $n$.  En el caso de $n = 2$: 
$$\Delta^2 f_i = \Delta(f_{i+1} - f_i) = f_{i+2} - 2f_{i+1} + f_i$$   
$$\nabla^2 f_i = \nabla(f_i - f_{i-1}) = f_i - 2f_{i-1} + f_{i-2}$$   
$$\delta^2 f_i = \delta\left(f_{i+\frac{1}{2}} - f_{i-\frac{1}{2}}\right) = f_{i+1} - 2f_i + f_{i-1}$$   
$$\Delta\nabla f_i = \Delta(f_i - f_{i-1}) = f_{i+1} - 2f_i + f_{i-1}$$   
$$\nabla\Delta f_i = \nabla(f_{i+1} - f_i) = f_{i+1} - 2f_i + f_{i-1}$$   

Conviene observar que en las ecuaciones anteriores las últimas tres diferencias son idénticas.  De hecho, podemos escribir la relación de identidad: $\delta^2 = \Delta\nabla = \nabla\Delta$. 

Las aproximaciones por diferencias se obtienen al aproximar los operadores diferenciales mediante los operadores de diferencias. 
$$\frac{d}{dx} \approx \frac{\Delta}{\Delta x}$$ 
$$\frac{d}{dx} \approx \frac{\nabla}{\nabla x}$$   
$$\frac{d}{dx} \approx \frac{\delta}{\delta x}$$ 

## Aproximación de derivadas parciales por diferencias 
Las aproximaciones por diferencias para el caso de las derivadas parciales de funciones multidimensionales son esencialmente iguales al caso de la diferenciación numérica de las funciones unidimensionales. 

Consideremos la función $f(x,y)$.  La aproximación por diferencia es $f_x = \frac{\partial f}{\partial x}$ en $x = x_0$ y $y = y_0$.  por lo tanto: 
$$f_x \approx \frac{f(x_0 + \Delta x, y_0) - f(x_0, y_0)}{\Delta x}$$   
$$f_x \approx \frac{f(x_0 + \Delta x, y_0) - f(x_0 - \Delta x, y_0)}{2\Delta x}$$   
$$f_x \approx \frac{f(x_0, y_0) - f(x_0 - \Delta x, y_0)}{\Delta x}$$   

Las aproximaciones por diferencias centrales para las segundas derivadas parciales de $f(x,y)$ en $x_0$ y $y_0$ son: 
$$f_{xx} \approx \frac{f(x_0 + \Delta x, y_0) - 2f(x_0, y_0) + f(x_0 - \Delta x, y_0)}{(\Delta x)^2}$$   
$$f_{yy} \approx \frac{f(x_0, y_0 + \Delta y) - 2f(x_0, y_0) + f(x_0, y_0 - \Delta y)}{(\Delta y)^2}$$   
$$f_{xy} \approx \frac{f(x_0 + \Delta x, y_0 + \Delta y) - f(x_0 + \Delta x, y_0 - \Delta y)}{4\Delta x\Delta y} + \frac{-f(x_0 - \Delta x, y_0 + \Delta y) + f(x_0 - \Delta x, y_0 - \Delta y)}{4\Delta x\Delta y}$$   


**Ejemplo 5.2**

La tabla siguiente muestra los valores de una función bidimensional $f(x, y)$:

| $y \backslash x$ | $1.0$ | $1.5$ | $2.0$ | $2.5$ | $3.0$ |
| --- | --- | --- | --- | --- | --- |
| **1.0** | 1.63 | 2.05 | 2.50 | 2.98 | 3.49 |
| **1.5** | 1.98 | 2.51 | 3.08 | 3.69 | 4.33 |
| **2.0** | 2.28 | 2.91 | 3.61 | 4.37 | 5.17 |
| **2.5** | 2.64 | 3.25 | 4.08 | 5.00 | 5.98 |
| **3.0** | 2.65 | 3.50 | 4.48 | 5.57 | 6.76 |

a) Utilice las aproximaciones por diferencias centrales para evaluar las siguientes derivadas parciales:

$$f_x(2, 2), \quad f_y(2, 2), \quad f_{yy}(2, 2) \quad \text{y} \quad f_{xy}(2, 2)$$

b) Use la aproximación por diferencias hacia adelante de tres puntos para evaluar las siguientes derivadas parciales:

$$f_x(2, 2), \quad f_y(2, 2)$$

**(Solución)**

Al emplear la definición $\Delta x = \Delta y = h = 0.5$, se hacen los cálculos como sigue:

a)

$$f_x(2, 2) = \frac{f(2 + h, 2) - f(2 - h, 2)}{2h} = \frac{4.37 - 2.91}{(2)(0.5)} = 1.46$$

$$f_y(2, 2) = \frac{f(2, 2 + h) - f(2, 2 - h)}{2h} = \frac{4.08 - 3.08}{(2)(0.5)} = 1.00$$

$$f_{yy}(2, 2) = \frac{f(2, 2 + h) - 2f(2, 2) + f(2, 2 - h)}{h^2} = \frac{4.08 - 2(3.61) + 3.08}{(0.5)^2} = -0.24$$

$$f_{xy}(2, 2) = \frac{f(2 + h, 2 + h) - f(2 + h, 2 - h) - f(2 - h, 2 + h) + f(2 - h, 2 - h)}{(2h)^2}$$

$$= \frac{5.0 - 3.69 - 3.25 + 2.51}{[2(0.5)]^2} = 0.57$$

b)

$$f_x(2, 2) = \frac{-f(2 + 2h, 2) + 4f(2 + h, 2) - 3f(2, 2)}{2h} = \frac{-(5.17) + 4(4.37) - 3(3.61)}{(2)(0.5)} = 1.48$$

$$f_y(2, 2) = \frac{-f(2, 2 + 2h) + 4f(2, 2 + h) - 3f(2, 2)}{2h} = \frac{-(4.48) + 4(4.08) - 3(3.61)}{(2)(0.5)} = 1.01$$