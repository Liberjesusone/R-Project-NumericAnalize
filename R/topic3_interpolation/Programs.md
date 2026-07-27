# PROGRAMA 2-1 Interpolación de Lagrange

**A) Explicaciones**

Este programa interpola una tabla de valores de funciones mediante una fórmula de interpolación de Lagrange. El usuario define la tabla de valores en los enunciados de los datos. Después de ejecutar la definición, la computadora pide el valor de *x* para el que se evaluará la fórmula de interpolación.

Aunque el programa está diseñado para la interpolación, se puede utilizar también para la extrapolación. Sin embargo, en este caso se imprime un mensaje “*X* está en el rango de extrapolación”.

**B) Variables**

**K:** número de datos
**F(I), X(I):** datos dados
**YRES:** valor numérico de la interpolación para un valor dado de *x*
**XA:** valor de *x* para el cual hay que evaluar la fórmula de interpolación

**C) Listado**

```fortran
C----CSL/F2-1.FOR      INTERPOLACION DE LAGRANGE
      DIMENSION F(0:10),X(0:10)
C           N ES EL ORDEN DEL POLINOMIO DE INTERPOLACION
      DATA N /3/
      DATA (X(I),I=0,3) /  1.,    2.,   3.,   4./
      DATA (F(I),I=0,3) / .671,.620,.567,.512 /
      PRINT *
      PRINT *,'CSL/F2-1      INTERPOLACION DE LAGRANGE'
      PRINT *
      PRINT *, 'TABLA DE VALORES UTILIZADOS'
      PRINT *, '---------------------------------------------------'
      PRINT *, '          I       X(I)            F(I)  '
      DO 37 I=0,N
        PRINT *, I,X(I),F(I)
37    CONTINUE
      PRINT *, '---------------------------------------------------'
45    PRINT *, 'DAR X ?'
      READ *,XA
      IF (XA.LT.X(0). OR . XA.GT.X(N)) PRINT *,
     &  '      ADVERTENCIA: X ESTA EN EL RANGO DE EXTRAPOLACION'
      YRES=0
      DO I=0,N
        Z=1.0
        DO J=0,N
          IF (I.NE.J) Z=Z*(XA-X(J))/(X(I)-X(J))
        END DO
        YRES = YRES + Z*F(I)
      END DO
      PRINT 200,XA, YRES
200   FORMAT('   RESULTADO DE LA INTERPOLACION :   G(',1PE12.5,') =',1PE12.5)
      PRINT *
      PRINT*,'   OPRIMA 1 PARA CONTINUAR, O 0 PARA TERMINAR
      READ *,K
      IF(K.EQ.1) GOTO 45
      PRINT*
      END
```

**D) Ejemplo de salida**

```
CSL/F2-1         INTERPOLACION DE LAGRANGE

TABLA DE VALORES UTILIZADOS
---------------------------------------------------
          I       X(I)            F(I)
          0      1.000000        0.6710000
          1      2.000000        0.6200000
          2      3.000000        0.5670000
          3      4.000000        0.5120000
---------------------------------------------------

DAR X ?
3.66
   RESULTADO DE LA INTERPOLACION :   G( 3.66000E+00) = 5.30924E-01

DAR X ?
4.5
      ADVERTENCIA: X ESTA EN EL RANGO DE EXTRAPOLACION
   RESULTADO DE LA INTERPOLACION :   G( 4.50000E+00) = 4.83750E-01

DAR X ?
0.1
      ADVERTENCIA: X ESTA EN EL RANGO DE EXTRAPOLACION
   RESULTADO DE LA INTERPOLACION :   G( 1.00000E-01) = 7.15190E-01
```

---

# PROGRAMA 2-3 Tabla de diferencias divididas

**A) Explicaciones**

Este programa desarrolla una tabla de diferencias divididas. Todos los datos de entrada están definidos en las instrucciones DATA.

Antes de ejecutar el programa, el usuario debe definir la tabla de valores en las instrucciones DATA. Los valores muestra de las instrucciones DATA en el programa listado a continuación son del ejemplo 2.8. Las diferencias divididas se calculan en dos ciclos, uno para *K* y otro para *I*.

**B) Variables**

**NI:** número de puntos dados en la tabla de valores
**J:** número máximo de diferencias para cada *k*
**K:** orden de una diferencia
**K(I):** valores *x* de los puntos
**F(K, I):** diferencia dividida de orden *k*: F(0, *I*) es el valor de la función para el punto I.

**C) Listado**

```fortran
C----CSL/F2-3.FOR      TABLA DE DIFERENCIAS DIVIDIDAS CON
C                      PUNTOS DE LA MALLA SEPARADOS DE MANERA NO UNIFORME
      DIMENSION F(0:10,0:10), X(0:10)
      PRINT *
      PRINT *,'CSL/F2-3      TABLA DE DIFERENCIAS DIVIDIDAS '
      DATA NI/6/
      DATA (X(I), I=0,6)/0.1, 0.2, 0.4, 0.7, 1.0, 1.2, 1.3/
      DATA (F(I,0),I=0,6)/.99750, .99002, .96040, .88120, .76520
     1      ,.67113, .62009/
      DO K=1,NI
        J=NI-K
        DO I= 0,J
          F(I,K)=(F(I+1,K-1)-F(I,K-1))/(X(I+K)-X(I))
        END DO
      END DO
      PRINT *
      PRINT *,' I    X(I)      F(I)      F(I,I+1)   F(I,I+2),..'
      DO I=0, NI
        J=NI-I
        PRINT 440,I,X(I), (F(I,K),K=0,J)
      END DO
      PRINT *
440   FORMAT (1X,I2,8F9.5)
      END
```

**D) Ejemplo de salida**

```
CSL/F2-3         TABLA DE DIFERENCIAS DIVIDIDAS
----------------------------------------------------------------------------------------------------
I    X(I)      F(I)      F(I,I+1), F(I,I+2),..
0   0.10000   0.99750 -0.07480 -0.24433   0.02089   0.01479 -0.00239   0.00128
1   0.20000   0.99002 -0.14810 -0.23180   0.03419   0.01215 -0.00085
2   0.40000   0.96040 -0.26400 -0.20444   0.04635   0.01122
3   0.70000   0.88120 -0.38667 -0.16737   0.05644
4   1.00000   0.76520 -0.47035 -0.13350
5   1.20000   0.67113 -0.51040
6   1.30000   0.62009
----------------------------------------------------------------------------------------------------
```