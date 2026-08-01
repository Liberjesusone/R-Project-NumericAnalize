## <font color="#00b0f0">PROGRAMA 3-2 Búsqueda de raíces</font>

**A) Explicaciones**
Este programa busca intervalos que contengan raíces de una función. La entrada consiste en los límites inferior y superior de $x$ para la búsqueda y un tamaño del intervalo $h$. Los intervalos en donde cambia el signo de la función (que contienen una o un número impar de raíces) se imprimen. El programa FUNC define la ecuación a resolver y se puede cambiar para los problemas propios del lector.

Antes de ejecutar el programa, el usuario debe definir la ecuación a resolver en FUNC, que por el momento tiene el problema del ejemplo 3.2. Cuando se ejecuta el programa, la computadora pide en forma interactiva tres parámetros de entrada, A, B y H. Así, el valor de X se designa como el límite inferior para la búsqueda y el contador de intervalos es I = 0. En el subprograma se calcula el valor de la función para el valor actual de $x$. Para I = 1 el programa brinca a S-40. Para I > 1, se verifica el producto de Y y YB, donde YB es el valor de Y para el valor anterior de X. Si el producto es negativo, se imprimen los valores de X — H y X como un intervalo que contiene un número impar de raíces. En la línea siguiente a S-40, se remplaza YB por el actual Y; X se incrementa en H a continuación el programa va a S-20. El programa se detiene si X excede a B, el límite máximo.

**B) Variables**
* **A:** límite inferior de $x$ para la búsqueda
* **B:** límite superior de $x$ para la búsqueda
* **H:** tamaño de los intervalos, $h$
* **Y:** valor de la función en $x$
* **YB:** valor de la función en $x - h$
* **I:** contador de intervalos

**C) Listado**

```fortran
C-----CSL/F3-2.FOR     BUSQUEDA DE RAICES
      print *
      PRINT *, 'CSL/F3-2    BUSQUEDA DE RAICES '
1     PRINT *
      PRINT *,'INICIAL X ? '
      READ *,A
      PRINT *,'FINAL X ? '
      READ *,B
      PRINT *,'INCREMENTO DE X ?'
      READ *,H
      PRINT *
      I=0                        ! Inicialización del número de intervalos
      X=A
20    I=I+1
      IF (X .GT. B) THEN
        GOTO 45
      ELSE
        Y = FUNC(X)
        IF (I.EQ.1.OR.Y*YB .GT. 0) GOTO 40
        PRINT *
        PRINT 90, X-H,X
      END IF
90    FORMAT(' UN INTERVALO QUE PUEDE CONTENER UNA RAIZ:  [',
     #       F10.6,',',F10.6,']')
40    YB=Y
      X=X+H
      GOTO 20
45    PRINT *
      PRINT *, '                *** FIN DE LA BUSQUEDA '
      PRINT *
      PRINT *,'OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR '
      READ *,KS
      IF(KS.EQ.1) GO TO 1
      STOP
      END
C*************************
      FUNCTION FUNC(X)           ! Define la ecuación a resolver.
      FUNC = -19*(X-.5)*(X-1)+EXP(X)-EXP(-2*X)
      RETURN
      END

```

**D) Ejemplo de salida**

```text
CSL/F3-2        BUSQUEDA DE RAICES

INICIAL X ? 
-10
FINAL X ? 
10
INCREMENTO DE X ? 
1
UN INTERVALO QUE PUEDE CONTENER UNA RAIZ:  [   0.000000,   1.000000]
UN INTERVALO QUE PUEDE CONTENER UNA RAIZ:  [   1.000000,   2.000000]
UN INTERVALO QUE PUEDE CONTENER UNA RAIZ:  [   6.000000,   7.000000]
                *** FIN DE LA BUSQUEDA 

INICIAL X ? 
0
FINAL X ? 
1
INCREMENTO DE X ? 
0.001
UN INTERVALO QUE PUEDE CONTENER UNA RAIZ:  [   0.405998,   0.406998]
                *** FIN DE LA BUSQUEDA 

```



## <font color="#00b0f0">PROGRAMA 3-5 Método de Newton</font>
**A) Explicaciones**
Este programa calcula una raíz real con una estimación inicial.
La ecuación a resolver y su primera derivada se definen en el subprograma FUNC. Se da como entrada la estimación inicial para $x$. En S-20, el contador de iteración se incrementa en uno. En cada iteración se encuentran X y XD llamando a FUNC, y a continuación se actualiza X mediante el método de Newton. La iteración termina si la diferencia entre dos valores consecutivos de $x$ es menor que la tolerancia especificada como entrada; el programa se detiene.

**B) Variables**

**X:** valor de $x$
**XB:** valor anterior de $x$
**Y:** valor de $y$ para el valor actual de $x$
**YD:** $y'$ para el valor actual de $x$
**I:** contador de pasos de iteración

**C) Listado**

```fortran
C-----CSL/F3-5.FOR          ESQUEMA DE NEWTON
C   LA ECUACION A RESOLVER Y SU DERIVADA SE DEFINEN
C   EN LA SUBRUTINA FUNC
      PRINT*
      PRINT*, 'CSL/F3-5    ESQUEMA DE NEWTON '
      PRINT*
      PRINT*, 'TOLERANCIA ?'
      READ *,EP
      PRINT *
5     PRINT*, 'ESTIMACION INICIAL PARA LA RAIZ ? '
      READ *, X
      XB=X
      I=0
      PRINT *
      PRINT *,' IT.NO. N   X(N-1)         Y(N-1)         X(N) '
20    I=I+1
      CALL FUNC(X,Y,YD)
      X = X - Y/YD               !-- Esquema de Newton: encuentra la nueva x.
      PRINT 30, I,XB,Y,X
30    FORMAT(1X,I5,3X,1P4E14.6)
      IF (ABS(X-XB).GE.EP) THEN  !-- Prueba de convergencia
        XB=X
        GO TO 20
      END IF
      PRINT *
40    PRINT*,'-----------------------------------------'
      PRINT*,'   SOLUCION FINAL=',X
      PRINT*,'-----------------------------------------'
      PRINT *
      PRINT*
      PRINT*,' PARA CONTINUAR, OPRIMA 1'
      READ *,K
      IF(K.EQ.1) GOTO 5
      PRINT*
      END
C*********************************
      SUBROUTINE FUNC(X,Y,YD)    !-- Calcula y y y'.
      Y=X**3 - 5.0*X**2 + 6.*X
      YD=3.0*X**2-10.0*X + 6.
      RETURN
      END

```

**D) Ejemplo de salida**

```text
CSL/F3-5        ESQUEMA DE NEWTON

TOLERANCIA?
0.00001
ESTIMACION INICIAL PARA LA RAIZ?
4.0
  IT.NO. N    X(N-1)          Y(N-1)          X(N)
      1     4.000000E+00    8.000000E+00    3.428571E+00
      2     3.428571E+00    2.099123E+00    3.127820E+00
      3     3.127820E+00    4.508972E-01    3.017077E+00
      4     3.017077E+00    5.240059E-02    3.000376E+00
      5     3.000376E+00    1.127243E-03    3.000000E+00
      6     3.000000E+00    1.907349E-06    3.000000E+00
-----------------------------------------
   SOLUCION FINAL =    3.000000
-----------------------------------------

ESTIMACION INICIAL PARA LA RAIZ ?
1.4
  IT.NO. N    X(N-1)          Y(N-1)          X(N)
      1     1.400000E+00    1.344000E+00    2.033962E+00
      2     2.033962E+00   -6.673241E-02    1.999361E+00
      3     1.999361E+00    1.277924E-03    2.000000E+00
      4     2.000000E+00    9.536743E-07    2.000000E+00
-----------------------------------------
   SOLUCION FINAL=     2.000000
-----------------------------------------

```