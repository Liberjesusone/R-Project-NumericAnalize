## Integracion ́

```
Carlos Echeverria
```
Universidad de ”Los Andes”.
Facultad de Ingenier ́ıa.

```
March 7, 2022
```

## Integracion num ́ erica ́

```
Regla del trapecio
Regla 1/3 de Simpson
Regla 3/8 de Simpson
Formula de Newton–Cotes ́
Cuadratura de Gauss
```

## Regla del trapecio

```
Se obtiene al integrar la formula de interpolaci ́ on lineal. Se ́
escribe en la forma siguiente:
```
### I=

```
∫b
```
```
a
```
```
f(x)dx=
b−a
2
[f(a) +f(b)] +E (1)
```

## Regla del trapecio

```
La ecuacion se puede extender a varios intervalos y se puede ́
aplicarNveces, con una separacion uniforme ́ h
```
### I=

```
∫b
```
```
a
```
```
f(x)dx=
b−a
2
```
### [

```
f(a) + 2
```
```
N− 1
```
# ∑

```
j= 1
```
```
f(a+jh) +f(b)
```
### ]

### +E (2)

```
dondeh= (b−a)/N
```

## Regla del trapecio


## Regla del trapecio


## Regla del trapecio


## Regla del trapecio

```
El error de la regla del trapecio se define como
```
### E=

```
∫b
```
```
a
```
```
f(x)dx−
b−a
2
```
```
[f(a) +f(b)] (3)
```
```
El desarrollo en serie de Taylor def(x),f(a),f(b)en torno a
x ̄= (a+b)/2, con la hipotesis de que ́ fes anal ́ıtica en
a 6 x 6 b
f(x) =f(x ̄) +zf′(x ̄) +
z^2
2
```
```
f′′( ̄x) +···
```
```
dondez=x−x ̄. El primer termino de la ecuaci ́ on ́
∫b
```
```
a
```
```
f(x)dx=
```
```
∫h/ 2
```
```
−h/ 2
```
### [

```
f(x ̄) +zf′(x ̄) +
```
```
z^2
2
f′′( ̄x) +···
```
### ]

```
dz (4)
```

## Regla del trapecio

```
Al integrar obtenemos:
∫b
```
```
a
```
```
f(x)dx=hf( ̄x) +
```
### 1

### 24

```
h^3 f′′( ̄x) +··· (5)
```
```
Por otro lado el segundo termino ́
```
```
b−a
2
```
```
[f(a) +f(b)] =
h
2
```
### [

```
f(x ̄)−
h
2
```
```
f′( ̄x) +
```
### 1

### 2

```
h^2
4
```
```
f′′( ̄x)−···
```
```
+ f( ̄x) +
```
```
h
2
f′(x ̄) +
```
### 1

### 2

```
h^2
4
f′′(x ̄) +···
```
### ]

```
= hf( ̄x) +
```
### 1

### 8

```
h^3 f′′(x ̄)−··· (6)
```
```
y se obtiene
```
### E=

```
∫b
```
```
a
```
```
f(x)dx−
```
```
b−a
2
[f(a) +f(b)]'−
```
### 1

### 12

```
h^3 f′′( ̄x) (7)
```

## Regla del trapecio

```
Supongamos que tenemosNintervalos en[a,b], entonces
```
### E'−

### 1

### 12

```
(b−a)^3
N^3
```
```
N
```
# ∑

```
i= 1
```
```
f′′(x ̄i) (8)
```
```
ya queh= (b−a)/Ny ̄f′′=N−^1
```
```
N
```
# ∑

```
i= 1
```
```
f′′( ̄xi). La ecuacion queda ́
```
### E'−

### 1

### 12

```
(b−a)h^2 ̄f′′(x ̄i) (9)
```
```
valida si la funci ́ on ́ f(x)es anal ́ıtica en el intervalo.
```

## Regla del trapecio

```
Supongamos dos resultadosIhyI 2 hy sus errores son
```
```
Eh'Ch^2 y E 2 h'C( 2 h)^2 =C 4 h^2. (10)
```
```
Por otra parte el valor exacto se puede escribir como
I=Ih+Eh=I 2 h+E 2 h, de lo cual obtenemos
```
```
Eh−E 2 h=I 2 h−Ih. (11)
```
```
Al despejarC
C=
```
### 1

### 3

```
h−^2 (Ih−I 2 h)
```
```
As ́ı la primera ecuacion es ́
```
```
Eh'
```
### 1

### 3

```
(Ih−I 2 h) (12)
```

## Regla del trapecio

```
De lo cual se obtiene una integral mas precisa:
```
```
I=Ih+Eh'Ih+
```
### 1

### 3

```
(Ih−I 2 h). (13)
```
```
No es exacta pero el error es del ordenh^4. Esta tecnica se ́
llama Integracion de Romberg ́
```

## Regla del trapecio


## Regla de 1/3 de Simpson

```
Se basa en la interpolacion polinomial cuadr ́ atica. El polinomio ́
de Newton hacia adelante ajustado a tres puntos,x 0 ,x 1 ,x 2.
```
### I=

```
∫b
```
```
a
```
```
f(x)dx=
h^3
3
```
```
[f(a) + 4 f(x ̄) +f(b)] +E (14)
```
```
dondeh= (b−a)/2 y ̄x= (a+b)/2.
```
```
E'−(h^5 / 90 )f′′′′( ̄x) (15)
```

## Regla de 1/3 de Simpson

```
La regla extendida es:
```
### I =

```
∫b
```
```
a
```
```
f(x)dx
```
### =

```
h^3
3
```
### 

```
f(a) + 4
```
```
N− 1
```
# ∑

```
(i impari=^1 )
```
```
f(a+ih) + 2
```
```
N− 2
```
# ∑

```
(i pari=^2 )
```
```
f(a+ih) +f(b)
```
### 

### +E

### (16)

```
dondeh= (b−a)/N. El error es
```
### E'−

### N

### 2

```
h^5
90
f′′′′(x ̄) = (a−b)
h^4
180
f′′′′( ̄x) (17)
```
```
Para un dominio fijo[a,b], el error es proporcional ah^4.
```

## Regla de 1/3 de Simpson


## Regla de 1/3 de Simpson

```
Al comparar los resultados anteriores con los del ejemplo 4.
se puede ver que la regla extendida de Simpson es mucho mas ́
precisa que la regla extendida del trapecio, utilizando el mismo
numero de intervalos. Por ejemplo, la exactitud de la regla ́
extendida del trapecio con 32 intervalos es equivalente a la de
la regla extendida de Simpson con tan solo 4 intervalos. El ́
error de la regla extendida de Simpson es proporcional ah^4 ,
por lo que es dosordenes m ́ as grande que el de la regla ́
extendida del trapecio.
```

## Regla de 1/3 de Simpson

```
Los desarrollos de Taylor paraf 0 yf 2 en torno dex 1 , o en forma
equivalente, ̄x= (a+b)/2, se escriben como
```
```
f 0 = f 1 −hf 1 ′+
```
### 1

### 2

```
h^2 f 1 ′′−
```
### 1

### 6

```
h^3 f 1 ′′′+
```
### 1

### 24

```
h^4 f 1 ′′′′−···
```
```
f 2 = f 1 +hf 1 ′+
```
### 1

### 2

```
h^2 f 1 ′′+
```
### 1

### 6

```
h^3 f 1 ′′′+
```
### 1

### 24

```
h^4 f 1 ′′′′+···
```
```
Sustituyendolo enI=
```
### 1

### 3

```
[f 0 + 4 f 1 +f 2 ] +Eobtenemos
```
```
I= 2 hf 1 +
```
### 1

### 3

```
h^3 f 1 ′′+
```
### 1

### 36

```
h^5 f 1 ′′′′+···+E. (18)
```
```
Por otro lado, el desarrollo de Taylor def(x)alrededor dex 1 es
```
```
f(x) = f 1 +zf 1 ′+
```
### 1

### 2

```
z^2 f 1 ′′+
```
### 1

### 6

```
z^3 f 1 ′′′+
```
### 1

### 24

```
z^4 f 1 ′′′′−···
```
```
dondex=x 1 +zoz=x−x 1
```

## Regla de 1/3 de Simpson

```
La integracion analitica de este desarrollo en ́ [a,b]da como
resultado
∫b
```
```
a
```
```
f(x)dx = 2 hf 1 +
```
### 1

### 3

```
h^3 f 1 ′′+
```
### 1

### 60

```
h^5 f 1 ′′′′+···
```
```
Restamos y truncamos despues del t ́ ermino principal, con lo ́
que el error de la ecuacion est ́ a dado aproximadamente por ́
```
### E'−

### 1

### 90

```
h^5 f 1 ′′′′ (19)
```
```
Una desventaja de la regla extendida de Simpson es que el
numero total de intervalos debe ser par. Por otro lado, la regla ́
de 3/8 de Simpson, se aplicaunicamente a un n ́ umero de ́
intervalos que sea multipo de tres. Por lo tanto, al combinar las ́
reglas de 1/3 y 3/8, se puede considerar el caso tanto par
como impar de intervalos.
```

## Regla de 3/8 de Simpson

```
Para un dominio[a,b]dividido en tres intervalos, se escribe
como
I=
```
```
∫b
```
```
a
```
```
f(x)dx=
```
### 3

### 8

```
h[f 0 + 3 f 1 + 3 f 2 +f 3 ] +E (20)
```
```
dondeh= (b−a)/3,fi=f(a+ih). El error
```
### E'−

### 3

### 80

```
h^5 f′′′′( ̄x) (21)
```
```
con ̄x= (a+b)/2.
```

## Formulas de Newton–Cotes ́

```
Los metodos de integraci ́ on num ́ erica que se obtienen al ́
integrar las formulas de interpolacion de Newton reciben el ́
nombre de formulas de Newton–Cotes. La regla del trapecio y
las dos reglas de Simpson son casos de las formulas de ́
Newton-Cotes, las cuales se dividen en formulas cerradas y ́
abiertas.
```
```
Escribimos las formulas cerradas de Newton-Cotes en la ́
forma:
```
```
∫b
```
```
a
```
```
f(x)dx=αh[w 0 f 0 +w 1 f 1 +w 2 f 2 +···+wNfN] +E (22)
```
```
dondeαywson las constantes que aparecen en la tabla y
```

## Regla del trapecio


## Formulas de Newton–Cotes ́

```
fn=f(xn), xn=a+nh h= (b−a)/n
```
```
Recibe el nombre de formula cerrada, debido a que el dominio ́
de integracion est ́ a cerrado por el primer y ́ ultimo dato. ́
```
```
La integracion se puede extender m ́ as all ́ a de los puntos ́
extremos (formulas abiertas). Dichas formulas se escriben ́
como
∫b
```
```
a
```
```
f(x)dx=αh[w 0 f 0 +w 1 f 1 +w 2 f 2 +···+wN+ 2 fN+ 2 ] +E (23)
```
```
dondeh= (b−a)/(N+ 2 ). Las constantesαywise listan en
la tabla, en dondew 0 ywN+ 2 se igualan a cero debido a que
corresponden a los extremos del dominio. Puesto quew 0 y
wN+ 2 se anulan,f 0 yfN+ 2 son datos ficticios, que en realidad
no son necesarios.
```

## Formulas de Newton–Cotes ́

```
Si comparamos una formula abierta con una cerrada utilizando, ́
para el mismo numero ́ Nde datos, el error de la formula ́
abierta es significativamente mayor que el de la formula ́
cerrada. Por otro lado, se pueden utilizar las formulas abiertas ́
cuando no se dispone de los valores de la funcion en los ́
l ́ımites de integracion. ́
```

## Formulas de Newton–Cotes ́


## Formulas de Newton–Cotes ́


## Formulas de Newton–Cotes ́


## Cuadratura de Gauss

```
Utilizan puntos de Legendre (raices de polinomios de
Legendre). Las cuadraturas de Gauss no se pueden utilizar
para integrar una funcion dada en forma de tabla con intervalos ́
de separacion uniforme debido a que los puntos de Legendre ́
no estan separados de esa manera; sin embargo, son m ́ as ́
adecuadas para integrar funciones analiticas. La ventaja de las
cuadraturas de Gauss es que su precision es mayor que la de ́
las formulas de Newton-Cotes.
```
```
El error de la regla del trapecio es proporcional af′′. Si se usa
esta regla para integrar cada una de las funciones
f= 1 ,x,x^2 ,x^3 ,...,entonces los resultados seran exactos para ́
f=1 yf=x, pero existiran errores para ́ x^2 y las potencias
superiores dex. El error de la regla de Simpson es
proporcional afiv, por lo que es exacta si integramos
f= 1 ,x,x^2 yx^3.
```

## Cuadratura de Gauss

```
Sin embargo, el siguiente ejemplo muestra que una integracion ́
numerica con dos puntos se puede hacer exacta para el caso ́
de los polinomios de orden tres, si se optimizan los valoresx
de los datos.
```

## Cuadratura de Gauss


## Cuadratura de Gauss

```
La formula general de la cuadratura de Gauss de orden ́ Nes
exacta cuando se integra un polinomio de orden menor o igual
a 2N−1.
```

## Cuadratura de Gauss

```
Las cuadraturas de Gauss difieren en forma significativa de las
formulas de Newton-Cotes ya que losNpuntos de la reticula
(llamados puntos de Gauss) se obtienen mediante las ra ́ıces
del polinomio de LegendrePN(x) =0, dondePN(x)es el
polinomio de Legendre de ordenN.
```
```
La cuadratura de Gauss que se extiende sobre el intervalo
[− 1 , 1 ]esta dada por ́
∫ 1
```
```
− 1
```
```
f(x)dx=
```
```
N
```
# ∑

```
k= 1
```
```
wkf(xk) (24)
```
```
dondeNes el numero de puntos de Gauss, Los ́ wison Los
pesos y lasxison los puntos de Gauss dados en la tabla.
```

## Cuadratura de Gauss


## Cuadratura de Gauss

```
Por ejemplo, siN=4, la ecuacion es ́
```
```
∫ 1
```
```
− 1
```
```
f(x)dx = 0. 34785 f(− 0. 86113 ) + 0. 65214 f(− 0. 33998 )
+ 0. 65214 f( 0. 33998 ) + 0. 34785 f( 0. 86113 )(25)
```
```
La formula de integraci ́ on de Gauss puede aplicarse a ́
cualquier intervalo arbitrario[a,b]con la transformacion ́
```
```
x=
2 z−a−b
b−a
```
### (26)

```
dondezes la coordenada original ena<z<byxes la
coordenada normalizada en− 16 x 6 1. La transformacion de ́
xenzes
z=
```
```
(b−a)x+a+b
2
```
### (27)


## Cuadratura de Gauss

```
Por medio de esta transformacion, la integral se puede escribir ́
como
```
```
∫b
```
```
a
```
```
f(z)dz=
```
```
∫ 1
```
```
− 1
```
```
f(z) (dz/dx)dx=
```
```
b−a
2
```
```
N
```
# ∑

```
k= 1
```
```
wkf(zk) (28)
```
```
dondedz/dx= (b−a)/2. Los valores dezkse obtienen al
sustituirxken la ecuacion (27) por los puntos de Gauss. ́
```
```
Por ejemplo, supongamos queN=2,a=0 yb=2. Puesto
que los puntos de GaussxkparaN=2 en la coordenada
normalizadax,− 16 x 6 1, son± 0 .57735 (de la tabla), los
puntos correspondientes enzson
```

## Cuadratura de Gauss

```
z 1 =
```
### 1

### 2

### [( 2 − 0 )(− 0. 57735 ) + 0 + 2 ] = 0. 42265 ]

```
z 2 =
```
### 1

### 2

### [( 2 − 0 )( 0. 57735 ) + 0 + 2 ] = 1. 57735 ]

```
La derivada esdz/dx= (b−a)/ 2 =1. Por lo tanto, la
cuadratura de Gauss se escribe como
∫b
```
```
a
```
```
f(z)dz=
```
```
∫ 1
```
```
− 1
```
```
f(z) (dz/dx)dx= [( 1 )f( 0. 42264 ) + ( 1 )f( 1. 57735 )]
(29)
```

