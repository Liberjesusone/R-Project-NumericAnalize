# Guía de Integración Numérica — Capítulo 4 (Nakamura)

> Universidad de Los Andes · CeSiMo · Análisis Numérico
> Programas 4-1 y 4-2 implementados en R

---

# PARTE I — FORMULARIO

## 1. Regla del Trapecio

**Simple** (un solo intervalo $[a,b]$):

$$I = \int_a^b f(x)\,dx \approx \frac{b-a}{2}\bigl[f(a) + f(b)\bigr]$$

**Extendida** ($N$ intervalos, $h = \dfrac{b-a}{N}$):

$$\boxed{\;I \approx \frac{h}{2}\left[f_0 + 2\sum_{j=1}^{N-1} f_j + f_N\right]\;}$$

Pesos: $\;1,\,2,\,2,\,2,\,\dots,\,2,\,1$

**Error:** $\;E \approx -\dfrac{1}{12}(b-a)\,h^2\,\bar{f}''\;\Rightarrow\; E = O(h^2)$

---

## 2. Regla $1/3$ de Simpson

**Simple** (dos intervalos, $h = \dfrac{b-a}{2}$):

$$I \approx \frac{h}{3}\bigl[f_0 + 4f_1 + f_2\bigr]$$

**Extendida** ($N$ **par**, $h = \dfrac{b-a}{N}$):

$$\boxed{\;I \approx \frac{h}{3}\left[f_0 + 4\!\!\sum_{i\ \text{impar}}\!\! f_i + 2\!\!\sum_{i\ \text{par}}\!\! f_i + f_N\right]\;}$$

Pesos: $\;1,\,4,\,2,\,4,\,2,\,\dots,\,4,\,1$

**Error:** $\;E \approx -\dfrac{(b-a)h^4}{180}f^{(iv)}\;\Rightarrow\; E = O(h^4)$

---

## 3. Regla $3/8$ de Simpson

**Simple** (tres intervalos, $h = \dfrac{b-a}{3}$):

$$\boxed{\;I \approx \frac{3h}{8}\bigl[f_0 + 3f_1 + 3f_2 + f_3\bigr]\;}$$

Pesos: $\;1,\,3,\,3,\,1$

**Error:** $\;E \approx -\dfrac{3}{80}h^5 f^{(iv)}$

> **Uso:** se aplica cuando $N$ es múltiplo de 3. En el Programa 4-1 se usa
> para "arrancar" cuando $N$ es impar (consume 3 intervalos y deja un número
> par para la regla $1/3$).

---

## 4. Integración de Romberg (extrapolación de Richardson)

Con dos resultados del trapecio, uno con paso $h$ y otro con paso $2h$:

$$E_h \approx \frac{1}{3}\bigl(I_h - I_{2h}\bigr)$$

$$\boxed{\;I \approx I_h + \frac{1}{3}\bigl(I_h - I_{2h}\bigr)\;}$$

**Efecto:** eleva la precisión de $O(h^2)$ a $O(h^4)$ **sin evaluar más puntos**.

---

## 5. Fórmulas cerradas de Newton-Cotes

$$\boxed{\;\int_a^b f(x)\,dx \approx \alpha\,h\left[w_0 f_0 + w_1 f_1 + \cdots + w_N f_N\right]\;}$$

donde $\;h = \dfrac{b-a}{N}$, $\;x_n = a + nh$, $\;\alpha = \dfrac{Q}{R}$.

### Tabla 4.2 — Constantes (fórmulas cerradas)

| $N$ | $\alpha$ | Pesos $w_i$ | Error |
|:---:|:---:|:---|:---|
| 1 | $1/2$ | $1,\ 1$ | $-\frac{1}{12}h^3 f''$ |
| 2 | $1/3$ | $1,\ 4,\ 1$ | $-\frac{1}{90}h^5 f^{(iv)}$ |
| 3 | $3/8$ | $1,\ 3,\ 3,\ 1$ | $-\frac{3}{80}h^5 f^{(iv)}$ |
| 4 | $2/45$ | $7,\ 32,\ 12,\ 32,\ 7$ | $-\frac{8}{945}h^7 f^{(vi)}$ |
| 5 | $5/288$ | $19,\ 75,\ 50,\ 50,\ 75,\ 19$ | $-\frac{275}{12096}h^7 f^{(vi)}$ |
| 6 | $1/140$ | $41,\ 216,\ 27,\ 272,\ 27,\ 216,\ 41$ | $-\frac{9}{1400}h^9 f^{(viii)}$ |
| 7 | $7/17280$ | $751,\ 3577,\ 1323,\ 2989,\ 2989,\ 1323,\ 3577,\ 751$ | $-\frac{8183}{518400}h^9 f^{(viii)}$ |
| 8 | $4/14175$ | $989,\ 5888,\ -928,\ \mathbf{10496},\ -4540,\ \mathbf{10496},\ -928,\ 5888,\ 989$ | $-\frac{2368}{467775}h^{11}f^{(x)}$ |
| 9 | $9/89600$ | $2857,\ 15741,\ 1080,\ 19344,\ \mathbf{5778},\ \mathbf{5778},\ 19344,\ 1080,\ 15741,\ 2857$ | $-\frac{173}{14620}h^{11}f^{(x)}$ |
| 10 | $5/299376$ | $16067,\ 106300,\ -48525,\ 272400,\ -260550,\ 427368,\ \dots$ (simétrico) | $-\frac{1346350}{326918592}h^{13}f^{(xii)}$ |

> ⚠️ **Erratas detectadas — tres en total.** Los valores correctos son los
> marcados en **negrita**:
> - $N=8$, peso: el libro imprime $10946$ → correcto **$10496$**
> - $N=8$, $\alpha$: la tabla teórica imprime $14/14175$ → correcto **$4/14175$**
> - $N=9$, peso: el libro imprime $5788$ → correcto **$5778$**
>
> La tercera es especial: los dos documentos del profesor **se contradicen**.
> `Integracion.md:282` (tabla) dice $14/14175$, pero `problems.md:186`
> (listado FORTRAN) dice `DATA .../4,14175/`. El código está bien; la tabla
> impresa está mal — por eso el programa original del libro funciona.
>
> Ver la sección "Erratas verificadas" al final para la demostración.

### Regla de verificación (¡memorizala!)

Toda fórmula de Newton-Cotes debe integrar exactamente $f(x)=1$. Eso obliga a:

$$\boxed{\;\alpha \cdot \sum_i w_i = N\;}$$

Si esto falla, hay un error en la tabla. Es la forma más rápida de detectar
un peso mal copiado.

---

## 6. Cuadratura de Gauss

$$\int_{-1}^{1} f(x)\,dx \approx \sum_{k=1}^{N} w_k\, f(x_k)$$

Para un intervalo arbitrario $[a,b]$, con el cambio de variable
$z = \dfrac{(b-a)x + a + b}{2}$:

$$\boxed{\;\int_a^b f(z)\,dz \approx \frac{b-a}{2}\sum_{k=1}^{N} w_k\, f(z_k)\;}$$

### Tabla 4.4 — Puntos y pesos de Gauss

| $N$ | $\pm x_i$ | $w_i$ |
|:---:|:---|:---|
| 2 | $0.577350269$ | $1.000000000$ |
| 3 | $0$ <br> $0.774596669$ | $0.888888889$ <br> $0.555555556$ |
| 4 | $0.339981043$ <br> $0.861136312$ | $0.652145155$ <br> $0.347854845$ |
| 5 | $0$ <br> $0.538469310$ <br> $0.906179846$ | $0.568888889$ <br> $0.478628670$ <br> $0.236926885$ |

**Propiedad clave:** la cuadratura de Gauss de orden $N$ es exacta para
polinomios de grado $\leq 2N-1$ — el doble de lo que logra Newton-Cotes
con el mismo número de puntos.

---

# PARTE II — TEORÍA

## ¿De dónde salen estas fórmulas?

Todas nacen de la misma idea: **reemplazar $f(x)$ por un polinomio
interpolante e integrar ese polinomio** (que sí sabemos integrar).

| Grado del polinomio | Regla que resulta | Puntos |
|---|---|---|
| Recta (grado 1) | Trapecio | 2 |
| Parábola (grado 2) | Simpson $1/3$ | 3 |
| Cúbica (grado 3) | Simpson $3/8$ | 4 |
| Grado $N$ | Newton-Cotes orden $N$ | $N+1$ |

Por eso el capítulo 4 depende directamente del capítulo 2 (interpolación):
**integrar numéricamente = interpolar y luego integrar exacto.**

---

## Deducción del error del trapecio

Se desarrolla $f(x)$ en serie de Taylor alrededor del punto medio
$\bar{x} = \frac{a+b}{2}$, con $z = x - \bar{x}$:

$$f(x) = f(\bar{x}) + z f'(\bar{x}) + \frac{z^2}{2}f''(\bar{x}) + \cdots$$

**Integrando el desarrollo** (valor verdadero):

$$\int_a^b f(x)\,dx = h\,f(\bar{x}) + \frac{1}{24}h^3 f''(\bar{x}) + \cdots$$

**Aplicando la regla del trapecio** al mismo desarrollo:

$$\frac{b-a}{2}[f(a)+f(b)] = h\,f(\bar{x}) + \frac{1}{8}h^3 f''(\bar{x}) + \cdots$$

**Restando:**

$$E = \left(\frac{1}{24} - \frac{1}{8}\right)h^3 f'' = -\frac{1}{12}h^3 f''(\bar{x})$$

Extendido a $N$ intervalos con $h = (b-a)/N$:

$$E \approx -\frac{1}{12}(b-a)\,h^2\,\bar{f}''$$

**Lo importante:** el error es proporcional a $h^2$. Si duplicás $N$
(mitad de $h$), el error se divide por **4**.

---

## ¿Por qué Simpson es tanto mejor?

El error de Simpson es $O(h^4)$: si duplicás $N$, el error se divide por **16**.

Comparación real (ejemplo 4.1, valor exacto $11.7286$):

| $N$ | Trapecio | Error | Simpson | Error |
|---:|---:|---:|---:|---:|
| 2 | 12.7627 | $1.03$ | 11.7809 | $-0.052$ |
| 4 | 11.9895 | $-0.26$ | 11.7318 | $-0.0032$ |
| 8 | 11.7940 | $-0.065$ | 11.7288 | $-0.0002$ |
| 16 | 11.7449 | $-0.016$ | 11.7286 | $0$ |
| 32 | 11.7326 | $-0.0040$ | 11.7286 | $0$ |

> **El dato que impresiona:** el trapecio con $N=32$ tiene la misma
> precisión que Simpson con $N=4$. Simpson necesita **8 veces menos
> evaluaciones** para la misma exactitud.

**¿Por qué gana Simpson?** Aunque la parábola es de grado 2, por simetría
la regla $1/3$ integra exactamente también los polinomios de grado 3 —
gana un orden "de regalo".

---

## Romberg: precisión gratis

Si sabés que $E_h \approx C h^2$, entonces con dos cálculos podés **estimar
y cancelar** el error:

$$E_h \approx C h^2, \qquad E_{2h} \approx C(2h)^2 = 4Ch^2$$

Como $I = I_h + E_h = I_{2h} + E_{2h}$:

$$E_h - E_{2h} = I_{2h} - I_h \;\Rightarrow\; -3Ch^2 = I_{2h} - I_h$$

$$\Rightarrow\; E_h \approx \frac{1}{3}(I_h - I_{2h})$$

$$\Rightarrow\; I \approx I_h + \frac{1}{3}(I_h - I_{2h})$$

**Verificado en el ejemplo 4.2:** con $I_{0.5}=11.9895$ e $I_{0.25}=11.7940$:

$$E_{0.25} = \tfrac13(11.7940 - 11.9895) = -0.0652$$
$$I = 11.7940 - 0.0652 = \mathbf{11.7288}$$

Ese resultado con solo 8 intervalos **iguala al trapecio con $N=128$**.

---

## Newton-Cotes: la generalización

Todas las reglas anteriores son casos particulares:

| Regla | $N$ | $\alpha$ | Pesos |
|---|:---:|:---:|---|
| Trapecio | 1 | $1/2$ | $1, 1$ |
| Simpson $1/3$ | 2 | $1/3$ | $1, 4, 1$ |
| Simpson $3/8$ | 3 | $3/8$ | $1, 3, 3, 1$ |

**Cerradas vs abiertas:**

- **Cerradas:** los extremos $a$ y $b$ **sí** se evalúan. Más precisas.
- **Abiertas:** los extremos **no** se evalúan ($w_0 = w_{N+2} = 0$). Sirven
  cuando $f$ es singular en los extremos (ej. $\int_0^1 \frac{\sin x}{x}dx$).

**Peligro de orden alto:** para $N \geq 8$ los pesos **alternan de signo** y
crecen mucho ($272400$, $-260550$, …). Sumar y restar números grandes produce
**cancelación catastrófica**. Por eso en la práctica no se usan órdenes altos;
se prefiere **subdividir** el intervalo y aplicar reglas de orden bajo
repetidamente (que es exactamente lo que hacen las reglas *extendidas*).

---

## Cuadratura de Gauss: la idea genial

Newton-Cotes **fija** los puntos equiespaciados y solo optimiza los pesos.
Gauss dice: *¿y si también elegimos dónde poner los puntos?*

Con $N$ puntos tenés $2N$ incógnitas ($N$ pesos + $N$ posiciones), así que
podés imponer $2N$ condiciones → exactitud hasta grado $2N-1$.

**Ejemplo 4.5 (deducción para $N=2$):** exigimos exactitud para
$f = 1, x, x^2, x^3$ en $[-1,1]$:

$$
\begin{aligned}
2 &= w_1 + w_2 \\
0 &= w_1 x_1 + w_2 x_2 \\
\tfrac23 &= w_1 x_1^2 + w_2 x_2^2 \\
0 &= w_1 x_1^3 + w_2 x_2^3
\end{aligned}
$$

Por simetría $x_2 = -x_1$, de donde $w_1 = w_2 = 1$ y:

$$x_1^2 = \tfrac13 \;\Rightarrow\; x_1 = \frac{1}{\sqrt3} = 0.577350269$$

Con **solo 2 puntos** integra exactamente cualquier cúbica. El trapecio
necesitaría infinitos.

**Limitación:** los puntos de Gauss no están equiespaciados, así que **no
sirve para datos tabulados** — solo para funciones que podés evaluar donde
quieras.

---

# PARTE III — EJEMPLOS RESUELTOS

## Ejemplo A — Trapecio a mano

Calcular $\int_0^2 \pi\left(1+\left(\frac{x}{2}\right)^2\right)^2 dx$ con $N=2$.

**Paso 1.** $h = \dfrac{2-0}{2} = 1$. Nodos: $x_0=0,\; x_1=1,\; x_2=2$.

**Paso 2.** Evaluar $f(x) = \pi(1+(x/2)^2)^2$:

| $i$ | $x_i$ | $(x/2)^2$ | $1+(x/2)^2$ | $f_i/\pi$ |
|:---:|:---:|:---:|:---:|:---:|
| 0 | 0 | 0 | 1 | $1$ |
| 1 | 1 | 0.25 | 1.25 | $1.5625$ |
| 2 | 2 | 1 | 2 | $4$ |

**Paso 3.** Aplicar la fórmula (pesos $1, 2, 1$):

$$I = \frac{h}{2}\bigl[f_0 + 2f_1 + f_2\bigr] = \frac{1}{2}\pi\bigl[1 + 2(1.5625) + 4\bigr]$$

$$I = 0.5\pi (8.125) = \mathbf{12.7627}$$

Error: $12.7627 - 11.7286 = 1.034$ (bastante grande con solo 2 intervalos).

---

## Ejemplo B — Simpson $1/3$ a mano

Mismo problema con $N=2$ (pesos $1, 4, 1$):

$$I = \frac{h}{3}\bigl[f_0 + 4f_1 + f_2\bigr] = \frac{1}{3}\pi\bigl[1 + 4(1.5625) + 4\bigr]$$

$$I = \frac{\pi}{3}(11.25) = \mathbf{11.7810}$$

Error: $0.052$ — **20 veces mejor** que el trapecio con los mismos 3 puntos.

Con $N=4$ ($h=0.5$, pesos $1,4,2,4,1$):

$$I = \frac{0.5}{3}\pi\bigl[1 + 4(1.0625)^2 + 2(1.25)^2 + 4(1.5625)^2 + 4\bigr] = \mathbf{11.7318}$$

---

## Ejemplo C — Simpson mixto ($N$ impar)

Integrar con $N=5$ en $[0,2]$. Como 5 es impar, no se puede usar $1/3$ sola.

**Estrategia del Programa 4-1:**

$$h = \frac{2-0}{5} = 0.4$$

1. **Regla $3/8$** sobre los 3 primeros intervalos $[0,\ 1.2]$:
   $$I_1 = \frac{3(0.4)}{8}\bigl[f_0 + 3f_1 + 3f_2 + f_3\bigr]$$
   con $x = 0,\ 0.4,\ 0.8,\ 1.2$

2. **Regla $1/3$** sobre los 2 restantes $[1.2,\ 2.0]$:
   $$I_2 = \frac{0.4}{3}\bigl[f_3 + 4f_4 + f_5\bigr]$$
   con $x = 1.2,\ 1.6,\ 2.0$

3. $I = I_1 + I_2 = \mathbf{11.73095}$ ✓ (coincide con el libro)

---

## Ejemplo D — Newton-Cotes orden 5

Calcular $\int_0^2 \sin(x)\,dx$ con 6 datos ($N = 5$).

**Paso 1.** $\alpha = \dfrac{5}{288}$, pesos $= 19, 75, 50, 50, 75, 19$,
$h = \dfrac{2}{5} = 0.4$.

**Paso 2.** Tabla:

| $j$ | $x_j$ | $f(x_j) = \sin x_j$ | $w_j$ | $w_j f_j$ |
|:---:|:---:|:---:|:---:|---:|
| 0 | 0.0 | 0.000000 | 19 | 0.0000 |
| 1 | 0.4 | 0.389418 | 75 | 29.2064 |
| 2 | 0.8 | 0.717356 | 50 | 35.8678 |
| 3 | 1.2 | 0.932039 | 50 | 46.6020 |
| 4 | 1.6 | 0.999574 | 75 | 74.9680 |
| 5 | 2.0 | 0.909297 | 19 | 17.2767 |
| | | | | $\sum = 203.9209$ |

**Paso 3.**

$$I = \alpha \cdot h \cdot \sum = \frac{5}{288}(0.4)(203.9209) = \mathbf{1.416117}$$

Exacto: $1 - \cos(2) = 1.416147$. Error: $3\times10^{-5}$ con solo 6 puntos.

---

# PARTE IV — CÓMO FUNCIONA CADA PROGRAMA

## Programa 4-1 (`program4_1.R`)

### Estructura

```
func()      → define f(x) a integrar
trapz()     → regla extendida del trapecio
simps()     → regla extendida de Simpson (maneja N par e impar)
romberg()   → extrapolación de Richardson
```

### `trapz()` — línea por línea

```r
trapz <- function(f, a, b, n) {
  h <- (b - a) / n              # espaciamiento uniforme
  s <- 0                        # acumulador
  for (i in 0:n) {              # recorre TODOS los nodos (n+1 puntos)
    x <- a + i * h              # posición del nodo i
    w <- if (i == 0 || i == n) 1 else 2   # ← peso: extremos 1, interiores 2
    s <- s + w * f(x)           # acumular w_i * f_i
  }
  s * h / 2                     # multiplicar por h/2 al final
}
```

**La clave:** el peso 2 en los interiores viene de que cada punto interior es
compartido por dos trapecios adyacentes — se cuenta dos veces.

### `simps()` — el caso difícil

```r
simps <- function(f, a, b, n) {
  h  <- (b - a) / n
  ss <- 0
  ls <- 0                       # ls = intervalos consumidos por la regla 3/8

  if (n %% 2 != 0) {            # ← N IMPAR: hay que arrancar con 3/8
    ls <- 3
    for (i in 0:3) {
      x <- a + h * i
      w <- if (i == 0 || i == 3) 1 else 3   # pesos 1,3,3,1
      ss <- ss + w * f(x)
    }
    ss <- ss * h * 3 / 8
    if (n == 3) return(ss)      # con N=3 ya terminamos
  }

  s <- 0
  m <- n - ls                   # intervalos que quedan (siempre PAR)
  for (i in 0:m) {
    x <- a + h * (i + ls)       # ← desplazado por ls
    w <- if (i %% 2 == 1) 4 else 2   # impar→4, par→2
    if (i == 0 || i == m) w <- 1     # extremos→1 (sobreescribe)
    s <- s + w * f(x)
  }
  ss + s * h / 3
}
```

**El truco de `ls`:** si $N$ es impar, se consumen 3 intervalos con la regla
$3/8$, y $N-3$ queda **par** — apto para la regla $1/3$. La variable `ls`
desplaza el índice para no volver a contar los nodos ya usados.

**El orden de los `if` importa:** primero se asigna 4 o 2 según paridad,
y **después** se sobreescribe con 1 si es extremo. Si fuera al revés, los
extremos quedarían mal.

### Salida verificada

```
Trapecio, N = 10 : 11.77048   (libro: 11.77047)
Simpson,  N =  5 : 11.73096   (libro: 11.73095)
```

> La diferencia en el último dígito es porque el FORTRAN usa `3.14159`
> y R usa `pi` con precisión completa.

---

## Programa 4-2 (`program4_2.R`)

### Estructura

```
NC_W, NC_Q, NC_R  → tabla 4.2 (pesos, numerador y denominador de alfa)
check_weights()   → verifica alfa*sum(w) = N para todos los órdenes
newton_cotes()    → aplica la fórmula general
demo_errata()     → demuestra las erratas de la tabla del libro
```

### Cómo se guarda la tabla

En el FORTRAN original, la matriz `W(0:20,10)` guarda **todo mezclado**:
los pesos en las posiciones $0..N$, y luego $Q$ y $R$ en las posiciones
$N+1$ y $N+2$:

```fortran
DATA (W(I,5),I=0,7)/19,75,50,50,75,19,  5,288/
                     └── pesos ──────┘  └Q┘└R┘
```

En R se separa en tres estructuras, que es mucho más legible:

```r
NC_W <- list("5" = c(19, 75, 50, 50, 75, 19))   # solo pesos
NC_Q <- c("5" = 5)                               # numerador de alfa
NC_R <- c("5" = 288)                             # denominador de alfa
```

### `newton_cotes()` — línea por línea

```r
newton_cotes <- function(f, a, b, k, verbose = TRUE) {
  n  <- k - 1                   # k = nº de DATOS, N = orden = k-1
  w  <- NC_W[[as.character(n)]] # pesos del orden N
  al <- NC_Q[ky] / NC_R[ky]     # alfa = Q/R
  h  <- (b - a) / n             # ← ojo: se divide por N, no por k

  ans <- 0
  for (j in 0:n) {
    x  <- a + j * h             # nodo j
    ans <- ans + f(x) * w[j+1]  # ← w[j+1] porque R indexa desde 1
  }
  ans * h * al                  # aplicar alfa*h al final
}
```

**Dos detalles que confunden:**

1. `k` es el **número de datos**, `n = k-1` es el **orden**. Con 6 datos el
   orden es 5.
2. R indexa vectores desde 1, pero los nodos van de $0$ a $N$. Por eso
   `w[j+1]` cuando el nodo es `j`.

### `check_weights()` — la validación

```r
# Toda fórmula de Newton-Cotes integra exactamente f(x) = 1:
#   I = alfa*h*sum(w) = alfa*((b-a)/N)*sum(w)  debe dar (b-a)
#   ⇒ alfa*sum(w) = N
```

Esta función corre al inicio y verifica los 10 órdenes. **Fue así como se
detectaron las erratas** de la tabla transcrita.

---

# PARTE V — ERRATAS VERIFICADAS

La tabla 4.2 transcrita del libro tiene **dos pesos mal impresos**.

## Detección

| $N$ | Valor del libro | $\sum w$ | $\alpha \cdot \sum w$ | ¿$= N$? |
|:---:|---|---:|---:|:---:|
| 8 | peso $10946$ | $29250$ | $8.254$ | ❌ |
| 8 | peso $\mathbf{10496}$ | $28350$ | $8.000$ | ✓ |
| 8 | $\alpha = 14/14175$ | $28350$ | $28.000$ | ❌ |
| 8 | $\alpha = \mathbf{4/14175}$ | $28350$ | $8.000$ | ✓ |
| 9 | peso $5788$ | $89620$ | $9.002$ | ❌ |
| 9 | peso $\mathbf{5778}$ | $89600$ | $9.000$ | ✓ |

### La errata del $\alpha$ en $N=8$

Esta es distinta a las otras dos: **los dos documentos del profesor se
contradicen entre sí**.

| Fuente | Dice | ¿Correcto? |
|---|---|:---:|
| `Integracion.md:282` (tabla teórica) | $14/14175$ | ❌ |
| `problems.md:186` (listado FORTRAN) | `DATA (W(I,8),I=9,10)/4,14175/` | ✅ |

Comprobación en el ejemplo 4.4 (longitud de arco, exacto $= 8$):

```
con  4/14175  ->  L = 8.000000    ✓
con 14/14175  ->  L = 28.000000   ✗  (factor 3.5 de error)
```

El código FORTRAN del libro está **bien**; lo que está mal es la tabla
impresa en el material teórico. Por eso el programa original funciona.

## Confirmación numérica

Usando el peso erróneo $5788$ en el ejemplo 4.4 (longitud de arco):

```
N=9 con peso 5778 (correcto): L = 8.00000   ← valor exacto
N=9 con peso 5788 (libro)   : L = 8.00198
              el libro reporta: L = 8.00197   ← coincide con el erróneo
```

## Consecuencia importante

El libro atribuye el crecimiento del error para $N \geq 9$ a **errores de
redondeo**. Eso **no es correcto**: con los pesos correctos, el método llega
a un error de $\sim10^{-13}$ en $N=10$ sin degradarse.

| $N$ | R (pesos correctos) | Libro |
|---:|---:|---:|
| 8 | 8.00000 | 8.00000 |
| 9 | 8.00000 | 8.00197 |
| 10 | 8.00000 | 7.99201 |

> El riesgo de **cancelación catastrófica** en Newton-Cotes de orden alto es
> real y sí existe — pero se manifiesta en órdenes mayores a los de la tabla,
> no en $N=9$.

---

# PARTE VI — GUION DE LA EXPOSICIÓN

## Lámina 1 — Portada
> "Capítulo 4 de Nakamura: Integración Numérica. Vamos a ver cómo calcular
> integrales definidas cuando no existe una antiderivada, implementando los
> programas 4-1 y 4-2 en R."

## Lámina 2 — Motivación
> "El problema es simple de enunciar: queremos $\int_a^b f(x)dx$. Pero muchas
> funciones **no tienen antiderivada elemental** — el ejemplo clásico es
> $e^{-x^2}$, que aparece en toda la estadística. Y muchas veces ni siquiera
> tenemos la función: tenemos una **tabla de datos** de un experimento.
> La idea central es: reemplazamos $f$ por un polinomio que sí sabemos
> integrar. Es decir, **interpolación del capítulo 2 + integración exacta**."

## Lámina 3 — Regla del Trapecio
> "La regla más simple: unimos los puntos con rectas y sumamos las áreas de
> los trapecios. Los pesos son 1 en los extremos y 2 en los interiores —
> el 2 aparece porque cada punto interior es compartido por dos trapecios.
> El error es $O(h^2)$: si duplico el número de intervalos, el error se
> divide por 4."

## Lámina 4 — Simpson
> "En vez de rectas usamos parábolas. La regla $1/3$ necesita $N$ par, y los
> pesos alternan $1, 4, 2, 4, 2, \ldots, 4, 1$. La $3/8$ usa cúbicas y necesita
> $N$ múltiplo de 3, con pesos $1, 3, 3, 1$. El error es $O(h^4)$ — al duplicar
> $N$ el error se divide por **16**, no por 4."

## Lámina 5 — Comparación de precisión
> "Acá está el resultado que más impresiona: el trapecio con 32 intervalos
> tiene la misma precisión que Simpson con 4. Simpson necesita **8 veces
> menos evaluaciones** de la función. Y eso importa mucho cuando cada
> evaluación es cara — por ejemplo una simulación completa."

## Lámina 6 — Análisis de error
> "Acá está el *por qué* de la lámina anterior. Cada regla tiene una fórmula
> de error, y todas dependen de una **derivada de orden alto** de $f$. El
> trapecio va con $f''$, Simpson con $f^{(iv)}$. Esa derivada mide cuánto se
> aleja la función real del polinomio que usamos — si $f$ **ya es** un
> polinomio de ese grado, la derivada vale cero y el error **desaparece**.
> Por eso Simpson integra exacto cualquier cúbica.
>
> Lo práctico es el exponente de $h$: con $O(h^2)$, duplicar $N$ divide el
> error por 4; con $O(h^4)$, lo divide por 16. Eso explica el resultado de
> la tabla anterior.
>
> Un detalle importante: la derivada se evalúa en un punto $\xi$ desconocido
> dentro del intervalo, así que estas fórmulas **acotan** el error, no lo dan
> exacto. Sirven para estimar cuántos intervalos necesito, no para corregir."

## Lámina 7 — Programa 4-1
> "Acá está la implementación. Lo interesante es cómo maneja el caso $N$
> impar: como la regla $1/3$ necesita $N$ par, el programa arranca aplicando
> $3/8$ sobre los primeros 3 intervalos, y lo que queda ya es par. La variable
> `ls` guarda cuántos intervalos consumió el arranque para no contarlos dos
> veces."

## Lámina 8 — Salida verificada
> "Las dos salidas del libro reproducidas en R: trapecio con $N=10$ da
> $11.77048$ contra $11.77047$ del libro, y Simpson con $N=5$ da $11.73096$
> contra $11.73095$. La diferencia en el último dígito es porque el FORTRAN
> usa la constante `3.14159` y R usa `pi` completo. El caso $N=5$ confirma que
> la combinación $3/8 + 1/3$ está bien implementada."

## Lámina 9 — Romberg
> "Esto es precisión gratis. Si sé que el error va como $h^2$, con dos
> cálculos puedo **estimar el error y restarlo**. Fíjense en el resultado:
> combinando los cálculos de $N=4$ y $N=8$ obtengo la misma precisión que el
> trapecio con $N=128$. Pasé de $O(h^2)$ a $O(h^4)$ sin evaluar ni un punto
> más."

## Lámina 10 — Newton-Cotes
> "Todas las reglas anteriores son casos particulares de una sola fórmula
> general. El trapecio es orden 1, Simpson $1/3$ es orden 2, Simpson $3/8$ es
> orden 3. La tabla 4.2 da los pesos hasta orden 10. Las fórmulas **cerradas**
> evalúan los extremos; las **abiertas** no, lo que sirve cuando la función es
> singular ahí."

## Lámina 11 — Programa 4-2 y erratas
> "Al implementarlo agregué una validación: toda fórmula de Newton-Cotes debe
> integrar exactamente $f=1$, lo que obliga a que $\alpha \cdot \sum w = N$.
> Corriendo esa verificación **detecté dos erratas en la tabla del libro**:
> en $N=8$ el peso $10946$ debería ser $10496$, y en $N=9$ el $5788$ debería
> ser $5778$. Y lo confirmé: usando el peso erróneo se reproduce exactamente
> el $8.00197$ que reporta el libro en el ejemplo 4.4. O sea que lo que el
> libro llama 'error de redondeo' es en realidad una errata de imprenta."

## Lámina 12 — Ejemplo 4.4 (longitud de arco)
> "Acá está Newton-Cotes aplicado a un problema real: la longitud de arco de
> una cardioide, cuyo valor exacto es exactamente 8. Con los pesos corregidos
> el método converge limpio hasta $10^{-13}$ en orden 10. Las dos filas en rojo
> de la columna del libro son justamente las que salen de los pesos con
> errata — es la confirmación práctica de lo de la lámina anterior."

## Lámina 13 — Cuadratura de Gauss
> "La idea genial: Newton-Cotes fija los puntos equiespaciados y solo optimiza
> los pesos. Gauss también optimiza **dónde poner los puntos**. Con $N$ puntos
> tengo $2N$ grados de libertad, así que logro exactitud hasta grado $2N-1$ —
> el doble. Con solo 2 puntos en $\pm 1/\sqrt3$ integra exactamente cualquier
> cúbica. La limitación: los puntos no son equiespaciados, así que no sirve
> para datos tabulados."

## Lámina 14 — Comparación de métodos
> "Este es el resumen del capítulo. Lo que quiero que quede: el trapecio es
> el más simple y robusto pero converge lento; Simpson es el caballo de
> batalla porque da $O(h^4)$ al mismo costo; Romberg exprime cálculos que ya
> hiciste; y Gauss da la máxima precisión pero solo sirve si podés evaluar la
> función donde quieras. Y la regla práctica de abajo es importante: **no
> subir el orden de Newton-Cotes**, mejor subdividir y repetir una regla
> simple — es el mismo argumento de los splines contra el polinomio único."

## Lámina 15 — Conclusiones
> "En resumen: trapecio para empezar o para datos ruidosos; Simpson es el
> caballo de batalla, mejor precisión por el mismo costo; Romberg cuando
> querés exprimir cálculos que ya hiciste; Gauss cuando podés evaluar la
> función donde quieras y necesitás máxima precisión. Y la lección
> transversal: **siempre validar las tablas** — la condición
> $\alpha\sum w = N$ me detectó dos erratas que llevaban décadas impresas."

---

# PARTE VII — RESUMEN COMPARATIVO

| Método | Orden error | Requiere | Ventaja | Limitación |
|---|:---:|---|---|---|
| **Trapecio** | $O(h^2)$ | Cualquier $N$ | Simplísimo, robusto | Lento a converger |
| **Simpson 1/3** | $O(h^4)$ | $N$ par | Mejor relación costo/precisión | $N$ debe ser par |
| **Simpson 3/8** | $O(h^4)$ | $N$ múltiplo de 3 | Completa a la $1/3$ | Poco usada sola |
| **Romberg** | $O(h^4)$ | Dos corridas del trapecio | Precisión "gratis" | Necesita $h$ y $2h$ |
| **Newton-Cotes** | hasta $O(h^{13})$ | Tabla de pesos | Generaliza todo | Cancelación en orden alto |
| **Gauss** | $2N-1$ exacto | Evaluar en puntos libres | Máxima precisión | No sirve para tablas |

## ¿Cuál usar?

```
¿Tengo la función o solo una tabla de datos?
├── Tabla con espaciado uniforme
│     └── Simpson 1/3 (o trapecio si N es impar y no querés complicarte)
├── Tabla con espaciado irregular
│     └── Trapecio (o splines + integración exacta)
└── Función evaluable en cualquier punto
      ├── Precisión moderada  → Simpson extendido
      ├── Máxima precisión    → Cuadratura de Gauss
      └── Ya tengo dos corridas → Romberg
```
