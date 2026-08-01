## Extrapolación vs interpolación
Imagina que tienes datos entre $x = 1$ y $x = 4$. Construyes $P(x)$ que pasa exactamente por esos 4 puntos.

- **Interpolación**: evalúas en $x = 3.66$ — está *dentro* del intervalo $[1,4]$. El polinomio tiene datos a ambos lados, "sabe" cómo comportarse. Resultado confiable.

- **Extrapolación**: evalúas en $x = 4.5$ o $x = 0.1$ — estás *fuera* del intervalo conocido. El polinomio no tiene ninguna restricción ahí, y los polinomios de grado alto tienden a dispararse hacia $\pm\infty$ fuera de sus nodos. Cuanto más te alejas, peor.

Por eso el programa imprime la advertencia pero sigue calculando: la extrapolación *puede* funcionar si te alejas poco, pero el riesgo crece rápido.

---

## Para qué sirve esto en la práctica

**Matemáticamente:**
- Las tablas de senos, cosenos y logaritmos que usaban los ingenieros antes de las calculadoras eran *exactamente* esto: se calculaban pocos valores exactos y se interpolaba entre ellos
- Integración numérica (trapecios, Simpson) — por dentro usan interpolación: aproximan $f(x)$ por un polinomio y lo integran de forma exacta
- Diferenciación numérica — misma idea al revés

**Computacionalmente:**
- **Gráficos**: las curvas suaves en programas de diseño (Illustrator, Blender, CAD) son splines — interpolación polinomial por trozos. Tú defines unos puntos de control y la computadora construye la curva que pasa por ellos
- **Animación**: un personaje tiene posición en el fotograma 1 y en el fotograma 30 — la computadora interpola los fotogramas intermedios
- **GPS y mapas**: los mapas digitales guardan puntos discretos del terreno y reconstruyen la superficie continua por interpolación
- **Señales de audio**: cuando subes audio digital a 44.1 kHz y lo reproduces, el chip tiene que reconstruir la onda continua — interpola entre las muestras
- **Física en videojuegos**: posiciones de objetos entre frames se interpolan

---

## Cómo lo hace la computadora

El problema en el fondo es este: dados $n+1$ puntos, existe **un único polinomio** de grado $\leq n$ que pasa por todos. La pregunta es cómo encontrarlo eficientemente.

**La forma obvia (y cara):** escribir el sistema lineal

$$\begin{pmatrix} 1 & x_0 & x_0^2 & \cdots & x_0^n \\ 1 & x_1 & x_1^2 & \cdots & x_1^n \\ \vdots & & & & \vdots \\ 1 & x_n & x_n^2 & \cdots & x_n^n \end{pmatrix} \begin{pmatrix} a_0 \\ a_1 \\ \vdots \\ a_n \end{pmatrix} = \begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix}$$

Eso se llama **matriz de Vandermonde**. Resolverlo cuesta $O(n^3)$ operaciones y es numéricamente inestable para $n$ grande. Nadie lo hace así.

**La forma de Lagrange (Prog 2-1):** construye directamente los polinomios "selectores" $L_k$ que valen $1$ en $x_k$ y $0$ en todos los demás nodos. No resuelve ningún sistema — directamente calcula el polinomio por construcción. Costo: $O(n^2)$.

**La forma de Newton con diferencias divididas (Prog 2-3):** construye el polinomio de forma *incremental*:

$$P(x) = \underbrace{f[x_0]}_{\text{orden 0}} + \underbrace{f[x_0,x_1]}_{\text{orden 1}}(x-x_0) + \underbrace{f[x_0,x_1,x_2]}_{\text{orden 2}}(x-x_0)(x-x_1) + \cdots$$

La clave: cada coeficiente $f[x_0,\ldots,x_k]$ se computa a partir del anterior con una resta y una división — la tabla triangular del Prog 2-3. Agregar un punto nuevo es solo añadir una columna. Costo: también $O(n^2)$, pero estructuralmente más flexible.

Los tres métodos (Vandermonde, Lagrange, Newton) dan el **mismo polinomio** — difieren solo en eficiencia y organización del cálculo. Eso es lo que hace que la interpolación sea tratable computacionalmente: encontrar representaciones inteligentes de un objeto matemático único.


# Presentación
Tu presentación tiene 9 slides. Te digo qué decir en cada una, qué es importante enfatizar y qué puede preguntar el profesor.

---

**Slide 0 · Título**
Solo di: *"El tema es Interpolación Polinomial, capítulo 2 de Nakamura. Vamos a ver qué es, los dos métodos principales del libro y sus implementaciones en R."* Nada más, pasa rápido.

---

**Slide 1 · ¿Qué es interpolación?**
Di esto: *"Tenemos n+1 puntos exactos de una función que quizás no conocemos en forma cerrada, o que es costosa de evaluar. La interpolación construye un polinomio que pasa exactamente por todos esos puntos, y lo usamos para estimar la función en puntos intermedios."*

Énfasis clave — el profe puede preguntar: **¿por qué no cualquier función, sino un polinomio?** Respuesta: los polinomios son fáciles de evaluar, derivar e integrar en la computadora. Son el bloque básico del análisis numérico.

También menciona la unicidad: *"El teorema garantiza que dado un conjunto de n+1 puntos con x distintos, hay exactamente un polinomio de grado ≤ n que los interpola."*

---

**Slide 2 · Fórmula de Lagrange**
Lee la fórmula pero explícala en palabras: *"$L_k(x)$ es un polinomio que vale exactamente 1 cuando $x = x_k$, y exactamente 0 cuando $x$ es cualquiera de los otros nodos. Entonces la suma $\sum L_k \cdot f_k$ selecciona el valor correcto en cada nodo y lo pondera suavemente entre ellos."*

Ese es el insight central — si el profe pregunta *"¿por qué funciona?"*, esa es la respuesta.

---

**Slide 3 · Implementación Prog 2-1**
Aquí describes el algoritmo en código: *"El Programa 2-1 implementa esto con dos lazos anidados. El lazo exterior recorre cada $i$, el interior calcula el producto de Lagrange $L_i$ excluyendo $j = i$. Al final acumula $Z \cdot f_i$."* Señala la tabla de datos del libro.

---

**Slide 4 · Resultados Lagrange**
Esto es importante para el profe: *"$x = 3.66$ está dentro del intervalo $[1,4]$ — interpolación legítima, el programa da $0.5309$. En cambio $x = 4.5$ y $x = 0.1$ están fuera — el programa lanza una advertencia pero sigue calculando. El resultado puede ser aceptable si nos alejamos poco, pero el error crece rápido."*

---

**Slide 5 · Diferencias divididas Prog 2-3**
Di: *"Newton reformula el polinomio en términos de diferencias divididas. La de orden 1 es simplemente la pendiente entre dos puntos. La de orden 2 es la 'aceleración'. Cada orden captura más curvatura de la función."*

La ventaja clave que SÍ va a preguntar: *"Si agrego un punto nuevo, con Lagrange tengo que recalcular todos los $L_k$ desde cero. Con Newton solo agrego una columna nueva a la tabla — los coeficientes anteriores no cambian."*

---

**Slide 6 · Error**
La fórmula del error es:
$$E(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \cdot \prod_{i=0}^{n}(x - x_i)$$

Explícala así: *"El error depende de dos cosas: qué tan 'curvada' es la función — eso es la derivada de orden n+1 — y qué tan lejos están los nodos del punto evaluado — eso es el producto. Por eso la extrapolación es peligrosa: el producto crece sin control fuera del intervalo."*

Menciona el **fenómeno de Runge**: *"Aumentar el grado del polinomio no siempre mejora la precisión. Para $\frac{1}{1+x^2}$ con nodos equiespaciados, aumentar $n$ empeora el error en los extremos. Por eso existen los nodos de Chebyshev."*

---

**Slide 7 · Comparación**
Resume en una frase: *"Los tres métodos producen el mismo polinomio único. Lagrange es conceptualmente más claro, Newton con diferencias divididas es más flexible para agregar puntos, y Newton hacia adelante/atrás es eficiente cuando los datos son equiespaciados."*

---

**Slide 8 · Resultados en R**
*"Implementé los dos programas del libro en R. El Prog 2-1 reproduce exactamente la salida del libro: $G(3.66) = 5.309 \times 10^{-1}$. El Prog 2-3 construye la tabla de diferencias divididas para $\cos(x)$ con 7 puntos no uniformes e interpola $\cos(0.5)$ con error del orden de $10^{-7}$."*

---

**Las dos preguntas más probables del profe:**

1. *"¿En qué se diferencia Lagrange de Newton si dan el mismo polinomio?"* → Organización del cálculo. Newton es incremental, Lagrange no.

2. *"¿Qué pasa si agrando n?"* → No necesariamente mejora. Fenómeno de Runge. Hay un n óptimo donde el error de truncamiento y el de redondeo se equilibran.


# Diferencia entre Newton hacia adelante y Newton con diferencias divididas
Son el mismo matemático, Newton, pero son **dos versiones distintas** de su método. La confusión viene de que los dos llevan su nombre.

---
## Newton hacia adelante (diferencias finitas Δ)
Requiere que los $x_i$ estén **equiespaciados** con paso $h$ constante. Define: $$\Delta f_i = f_{i+1} - f_i \qquad \Delta^2 f_i = \Delta f_{i+1} - \Delta f_i \qquad \cdots$$
El polinomio queda: $$P(x) = f_0 + s\,\Delta f_0 + \frac{s(s-1)}{2!}\Delta^2 f_0 + \cdots \qquad s = \frac{x - x_0}{h}$$
**Limitación**: solo funciona con datos de espaciado uniforme.

---
## Newton con diferencias divididas (Prog 2-3)
Funciona con **cualquier espaciado**. En vez de restar simplemente, *divide* por la separación entre nodos: $$f[x_i, x_j] = \frac{f[x_j] - f[x_i]}{x_j - x_i}$$
El polinomio queda: $$P(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \cdots$$

---
## Relación entre los dos
Son **el mismo polinomio** expresado distinto. Cuando el espaciado es uniforme, las diferencias divididas se reducen a diferencias finitas divididas por potencias de $h$: $$f[x_i, x_{i+1}] = \frac{\Delta f_i}{h} \qquad f[x_i, x_{i+1}, x_{i+2}] = \frac{\Delta^2 f_i}{2h^2}$$
O sea: Newton hacia adelante es un caso especial de diferencias divididas cuando $h$ es constante.

---
## ¿Cuál usa el Programa 2-3?
**Diferencias divididas** — mira la línea clave del código:
```r
# program2_3.R : 12
dd[i, k] <- (dd[i+1, k-1] - dd[i,k-1]) / (x_data[i+k-1] - x_data[i])
```

Ese denominador `x_data[i+k-1] - x_data[i]` es exactamente la división por la separación entre nodos — eso es la diferencia dividida. Por eso funciona con los datos de $\cos(x)$ que tienen espaciado irregular: $h = 0.1, 0.2, 0.3, 0.3, 0.2, 0.1$.

Si los datos fueran equiespaciados, ese denominador sería siempre una potencia de $h$ y coincidiría exactamente con Newton hacia adelante.


# Explicación del algoritmo 2-3

## Por qué cos(x) en el Prog 2-3
Mira el FORTRAN en `Programs.md` — los datos del libro **ya son** cos(x), solo que guardados como tabla:

```fortran
DATA (X(I), I=0,6) / 0.1, 0.2, 0.4, 0.7, 1.0, 1.2, 1.3 /
DATA (F(I,0),I=0,6) / .99750, .99002, .96040, .88120, .76520, .67113, .62009 /
```

Verifica: `cos(0.1) = 0.99750`, `cos(0.2) = 0.99002`... son exactamente esos valores. El libro los llama del "Ejemplo 2.8" que usa cos(x). El programa FORTRAN solo ve la tabla — no sabe que es cos. En R pusimos `cos(x)` explícitamente para poder comparar el resultado interpolado contra el valor exacto y calcular el error. El algoritmo es el mismo.

---

## Cómo funciona el algoritmo de diferencias divididas

El objetivo es construir esta tabla triangular — cada columna es un orden mayor:

```
I   x_i     f[·]       f[·,·]    f[·,·,·]   ...
0   0.1   0.99750   -0.07480   -0.24433   ...
1   0.2   0.99002   -0.14810   -0.23180   ...
2   0.4   0.96040   -0.26400   -0.20444   ...
3   0.7   0.88120   -0.38667   ...
4   1.0   0.76520   -0.47035
5   1.2   0.67113   -0.51040
6   1.3   0.62009
```

**Columna 1** — son los datos directos: $f[x_i] = f(x_i)$

**Columna 2** — pendiente entre pares adyacentes:
$$f[x_0, x_1] = \frac{f[x_1] - f[x_0]}{x_1 - x_0} = \frac{0.99002 - 0.99750}{0.2 - 0.1} = -0.07480$$

**Columna 3** — "pendiente de pendientes":
$$f[x_0, x_1, x_2] = \frac{f[x_1,x_2] - f[x_0,x_1]}{x_2 - x_0} = \frac{-0.14810 - (-0.07480)}{0.4 - 0.1} = -0.24433$$

Y así sucesivamente. El código hace exactamente eso:

```r
# program2_3.R : 12
dd[i, k] <- (dd[i+1, k-1] - dd[i,k-1]) / (x_data[i+k-1] - x_data[i])
```

`dd[i+1,k-1] - dd[i,k-1]` es la resta entre dos diferencias del orden anterior, y se divide por la separación entre los nodos extremos.

**Para evaluar** $P(x)$ en un punto nuevo, se usa la primera fila de la tabla como coeficientes:

$$P(x) = \underbrace{0.99750}_{f[x_0]} + \underbrace{(-0.07480)}_{f[x_0,x_1]}(x-0.1) + \underbrace{(-0.24433)}_{f[x_0,x_1,x_2]}(x-0.1)(x-0.2) + \cdots$$

Por eso en el código de evaluación el lazo siempre empieza desde `dd[1, ...]` — la primera fila es la que tiene todos los coeficientes.

La intuición: diferencia dividida de orden 1 = pendiente (como una derivada discreta). Orden 2 = curvatura. Orden 3 = cambio en curvatura. Cada columna captura un nivel más de variación de la función.

---
`term` es una productoria acumulada. Te explico por qué es necesaria.

Una vez tienes la primera fila de la tabla, el polinomio es:

$$P(x) = c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1) + c_3(x-x_0)(x-x_1)(x-x_2) + \cdots$$

donde $c_k = f[x_0, x_1, \ldots, x_k]$ son los valores de la primera fila.

El problema: cada término tiene un producto que crece — el término 3 necesita $(x-x_0)(x-x_1)(x-x_2)$, que incluye todo lo del término anterior más un factor nuevo.

En vez de recalcular ese producto desde cero en cada iteración, `term` lo va *acumulando*:

```r
result <- dd[1, 1]          # c₀  →  arranca con f[x₀]
term   <- 1                 # el producto empieza vacío (= 1)

# k=1:  term = 1 · (x - x₀)              result += term · c₁
# k=2:  term = (x-x₀) · (x - x₁)        result += term · c₂
# k=3:  term = (x-x₀)(x-x₁) · (x - x₂) result += term · c₃
for (k in 1:(n-1)) {
  term   <- term * (x_target - x_data[k])   # multiplica un factor más
  result <- result + term * dd[1, k+1]       # suma el término nuevo
}
```

Cada vuelta del lazo hace dos cosas: extiende el producto en un factor `(x - x_k)` y suma ese producto multiplicado por el coeficiente correspondiente de la primera fila.

Es la misma idea que sumar una serie donde cada término se obtiene del anterior multiplicando algo — en vez de calcular $n!$ desde cero en cada paso, llevas un acumulador que multiplica por $n$ en cada vuelta.


Porque el polinomio de Newton siempre arranca desde $x_0$.

Mira la estructura de la tabla:

```
fila 1:  dd[1,1]  dd[1,2]  dd[1,3]  dd[1,4]  ...
fila 2:  dd[2,1]  dd[2,2]  dd[2,3]  dd[2,4]  ...
fila 3:  dd[3,1]  dd[3,2]  dd[3,3]  ...
```

Lo que hay en cada posición es:

```
dd[1,1] = f[x₀]
dd[1,2] = f[x₀, x₁]
dd[1,3] = f[x₀, x₁, x₂]
dd[1,4] = f[x₀, x₁, x₂, x₃]

dd[2,1] = f[x₁]
dd[2,2] = f[x₁, x₂]
dd[2,3] = f[x₁, x₂, x₃]
```

El polinomio que se construye es:

$$P(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \cdots$$

Todos los coeficientes empiezan desde $x_0$ — por eso todos están en la fila 1. La columna sube porque sube el **orden** de la diferencia dividida, pero el punto de partida siempre es $x_0$.

Si usaras `dd[2, k+1]` estarías tomando coeficientes que arrancan desde $x_1$, que daría un polinomio diferente — válido matemáticamente, pero no el de Newton estándar que empieza en $x_0$.

Las filas 2, 3, 4... de la tabla son intermedias — se calculan para poder llegar a la fila 1, pero no se usan directamente en la evaluación final.