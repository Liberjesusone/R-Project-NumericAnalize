## Errores.

```
Carlos Echeverria
```
Universidad de ”Los Andes”.
Facultad de Ingenier ́ıa.

```
January 26, 2022
```

## Reglas del juego.

(^1) Son tres parciales, con un examen final.
(^2) La bibliograf ́ıa a emplear son:
Shoichiro Nakamura, Metodos num ́ ericos aplicados con ́
software. 1992.
https://www.dropbox.com/s/yyg4cauqb659z2y/Nakamura.pdf?dl=
S. C. Chapra y R. P. Canale. Metodos num ́ ericos para ́
ingenieros. Quinta Edicion. 2007. ́
https://www.dropbox.com/s/bvvqigyn05a0ljs/Chapra.pdf?dl=
Notas de clases
https://www.dropbox.com/s/n46da3v2gww1evg/Errores.pdf?dl=


## Introduccion a los c ́ alculos num ́ ericos. ́

```
El eficiente y preciso uso de los resultados numericos requiere ́
considerable experiencia y avanzados conocimientos para
evitar desatinos.
```
```
Una realista y exitosa solucion a un problema de ingenier ́ ́ıa
empieza con:
```
(^1) Un preciso modelo f ́ısico del problema y apropiada
comprension de las suposiciones empleadas. ́
(^2) Luego, transformarlo a un problema matematico y ́
usualmente se obtiene a traves de m ́ etodos num ́ ericos ́
que por definicion son aproximados ́
Seguidamente, esto determinara los objetivos, los datos ́
adecuados y sus verificaciones.


## Introduccion a los c ́ alculos num ́ ericos. ́

```
Una vez formulado el problema:
Deben disenarse los m ̃ etodos num ́ ericos, junto con un ́
analisis preliminar del error. ́
Un metodo num ́ erico que puede usarse para resolver el ́
problema se llama algoritmo. (Es un conjunto de
procedimientos y carente de ambiguedad). ̈
Ya decidido el algoritmo o conjunto de algoritmos, el
analista debe considerar todas las fuentes de error
(numero de iteraciones, vereficar la exactitud y dejar ́
espacio para la accion correctiva). ́
El programador debe transformar el algoritmo en un
conjunto de instrucciones detallado y sin ambiguedad. ̈
```

## Fuentes basicas de los errores. ́

```
Todas las soluciones numericas contiene errores de diferentes ́
or ́ıgenes, que las podemos resumir en:
Error del modelo o error del problema Es cuando las
hipotesis son muy alejadas de la realidad o sencillamente ́
estan erradas. Tambien, cuando las condiciones son ́
erradas.
Error del metodo ́ Cuando un problema en forma precisa
no puede resolverse en forma exacta o es muy dificil hallar
la solucion, se f ́ ormula una aproximaci ́ on del modelo. ́
Error residual Son los originados por las series infinitas,
al considerar solo una parte finita.
Error inicial Son los originados por los parametros cuyos ́
valores son conocidos aproximadamente
```

## Fuentes basicas de los errores. ́

```
Error de redondeo Originados por la representacion finita ́
de los numeros. ́
Error de operacion ́ Podemos distinguir dos casos:
operacion exacta de n ́ umeros aproximados ́
π− 1 / 6 ≈ 3. 1416 − 0. 1667 = 3. 2083
operacion aproximada de numeros exactos ́
0. 23 × 108 − 1. 36 × 10 −^8 ≈ 0. 23 × 108 , el valor verdadero es
22999999 .9999999864.
Error sistematico ́ Son aquellos que sin variar las
condiciones del ensayo, entra de igual modo en cada
resultado, puede ser originado por:
Defecto del instrumento
Las condiciones del ambiente
La metodolog ́ıa de la medicion ́
Precision limitada del instrumento ́
Las particularidades del experimento
```

## Fuentes basicas de los errores. ́

```
Error Causal o Accidental (fortuito) Son los que estan ́
vinculados con los factores que sufren pequenas ̃
variaciones (aleatorias) durante el experimento.
```
## Error absoluto, relativo y porcentual

```
Consideremos ”A” el valor exacto y ”a” el valor conocido que se
llamara aproximaci ́ on de ”A”. Definiciones: ́
Llamemos error absolutoξa=|A−a|
Llamenos error relativoδa=ξa/|A|
Llamenos error relativo porcentualεa=δa× 100 %
```

## Precision y exactitud ́


## Analisis de error en los m ́ etodos num ́ ericos ́

```
Existen dos causas principales de errores en los calculos ́
numericos: ́
```
(^1) **El error de truncamiento** : Se debe a las aproximaciones
utilizadas en la formula matematica del modelo. ́
(^2) **El error de redondeo** : Se asocia con el numero limitado ́
de d ́ıgitos con que se representa los numeros en una ́
computadora
**Series de Taylor**
f(x) = f(a)+h f′(a)+h
2
2 f
′′(a)+h^3
6 +f
′′′(a)+h^4
24 f
′′′′(a)
+h
5
5 !f
′′′′′(a)+···+hm
m!f
(m)(a)+··· (1)
dondeh=x−a


## Analisis de error en los m ́ etodos num ́ ericos ́

```
En las aplicaciones practicas, hay que truncar la serie ́
```
```
f(x) = f(a)+h f′(a)+h
2
2 f
```
```
′′(a)+h^3
6 +f
```
```
′′′(a)+h^4
24 f
```
```
′′′′(a)
```
```
+h
5
5 !f
```
```
′′′′′(a)+···+hN
N!f
```
(N)(a)+O(hN+ (^1) ) (2)
dondeO(hN+^1 )representa el error provocado por el
truncamiento de los terminos de orden ́ N+1 y superiores
O(hN+^1 ) =
hN+^1
(N+ 1 )!
f(N+^1 )(a+ξh), 0 ≤ξ≤ 1
o ́
O(hN+^1 ) =
hN+^1
(N+ 1 )!
f(N+^1 )(a) (3)


## Analisis de error en los m ́ etodos num ́ ericos ́

```
Bases de los n ́umeros
El extendido uso del sistema decimal oculta la existencia de
otros sistemas numericos. Sin embargo, en computaci ́ on se ́
usa otros sistemas numericos ́
Binario: base 2
Octal: base 8
Hexadecimal: base 16
Estos sistemas pueden traducirse con facilidad de uno al otro.
El octal y hexadecimal son mas cortas que en binario. El ́
hexadecimal proporciona un uso mas eficiente de la memoria. ́
La forma de expresarla es:( 3. 224 ) 10 ,( 1001. 11 ) 2 y
( 18 C 7. 90 ) 16
```
```
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
0 1 2 3 4 5 6 7 8 9 A B C D E F
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
El valor decimal de un numero en base r ́
```
```
(abcderfg.hijk)r
```
```
y se calcula
```
```
ar^6 +br^5 +cr^4 +dr^3 +er^2 +fr^1 +g+hr−^1 +ir−^2 +jr−^3 +kr−^4 (4)
```
```
N ́umero dentro del Hardware de la computadora
bit: del inglesBinary digit, representa un d ́ıgito binario.
Ademas, el bit es la unidad m ́ ́ınima de informacion empleada ́
en la teor ́ıa de la informacion. ́
byte u octeto: Generalmente se refiere a 8 bits. Formalmente
es una secuencia de bits contiguos, cuyo tama ̈ no depende del ̃
codigo de informaci ́ on o de car ́ acteres en que se est ́ e usando. ́
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Enteros: En el sistema de enumeracion binario ́
```
```
±akak− 1 ak− 2 ···a 2 a 1 a 0
y su valor decimal es
```
```
I=±[ak 2 k+ak− 12 k−^1 +···a 222 +a 12 +a 0 ]
Por ejemplo
± 110101
I= ±[ 1 x 25 + 1 x 24 + 0 x 23 + 1 x 22 + 0 x 2 + 1 ] =± 53
El valor dekse limita por el diseno del hardware. El primer bit ̃
registra el signo, 0→+y 1→−. Los restantes se usan para
losak
```
### 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

### 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15


## Analisis de error en los m ́ etodos num ́ ericos ́

```
El valor decimal de lo anterior es:
14
```
# ∑

```
i= 0
```
```
2 i= 32767
```
```
N ́umeros reales: Cuando se introduce un numero decimal, ́
primero se convierte al binario mas cercano ́
```
```
(± 0 .abbbbb...bbbb) 2 x 2 z (5)
dondeasiempre es 1, cadabes un d ́ıgito binario yzes un
exponente que se expresa en binario
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Los 32 bits se distribuyen de la manera siguiente
```
```
11111111 11111111 11111111 11111111
seeeeeee emmmmm mmmmm mmmmm
El exponente va desde 0 hasta 2^8 − 1 =250. Para registrar
exponentes positivos y negativos el exponente en decimal es
sesgado por 128 y despues convertido binario. ́
Por ejemplo si el exponente es−3, entonces− 3 + 128 =125,
se convierte en binario y se almacena en los 8 bits. Por lo
tanto, los exponentes que se pueden almecenar van desde
0 − 128 =−128 hasta 255− 128 = 127
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
En el caso hexadecimal
```
```
(± 0 .abbbbb) 16 x 16 k (6)
```
```
dondeaes un d ́ıgito hexadecimal distinto a cero,bes un d ́ıgito
hexadecimal ykes un esponente expresado en binario
```
```
11111111 11111111 11111111 11111111
seeeeeee mmmmm mmmmm mmmmm
La maxima mantisa positiva en hexadecimal ́
( 0 .FFFFFF) 16 = 1 −( 16 )−^6.
El maximo exponente 2 ́^7 − 1 − 64 =63.
Por lo tanto el maximo valor( 1 − 16 −^6 ) 1663 ' 7. 23 x 1075 ,
y el menor( 0. 100000 ) 16 x 16 −^64 = 5. 39 x 10 −^79
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Errores de redondeo al almacenar un n ́umero en memoria
```
```
Elepsilon de la m ́ aquina, ́ ε, es el tamano del intervalo entre 1 y ̃
el siguiente numero mayor que 1 distinguible de 1 ́
Se puede determinar con:
>epsilon <- 1.
>while(epsilon + 1 > 1) epsilon <-
epsilon/2.
>print(epsilon)
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Causas de efecto por redondeo
Considere la suma 1+ 0 .00001:
```
```
( 1 ) 10 = ( 0 .1000 0000 0000 0000 0000 0000) 2 x 21
( 0. 00001 ) 10 = ( 0 .1010 0111 1100 0101 1010 1100) 2 x 2 −^16
```
```
la suma de estos
```
### ( 1 ) 10 + ( 0. 00001 ) 10

### = ( 0 .1000 0000 0000 0000 0101 0011

```
1110 0010 1101 0110 0) 2 x 21
' ( 0 .1000 0000 0000 0000 0101 0100) 2 x 21 = 1. 0000100136
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Errores por redondeo implicados al restar numeros ́
```
```
lim
θ→ 0
```
```
f(x+θ)−f(x)
θ
```
```
=f′(x) (7)
```
```
Hacemosf(x) =sen(x)y calculamos
```
```
d=
```
```
sin( 1 +θ)−sin( 1 )
θ
```
### (8)

```
>for(k in 1:20) {
>... t = 10.ˆ(-k)
>... d = (sin(1.+t)-sin(1))/t
>... error = 0.54030 - d
>... print(k,t,d,error)}
```

## Analisis de error en los m ́ etodos num ́ ericos ́

```
Para analizar el redondeo en la resta consideremos
1. 00001 − 1
```
```
1. 00001 = ( 0 .1000 0000 0000 0000 0101 0100) 2 x 21
entonces,( 1. 00001 − 1 )es
( 0 .1000 0000 0000 0000 0101 0100) 2 x 21 −( 0. 1 ) 2 x 21
= ( 0 .0000 0000 0000 0000 0101 0100) 2 x 21
= ( 0 .1010 1) 2 x 2 −^16 = 1. 00136 x 10 −^5
```
### (9)

```
Al comparar con el valor exacto
( 0. 0000100136 − 0. 00001 )/ 0. 00001 = 0. 00136 o 0 ́. 136 %
```

