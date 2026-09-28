### PROGRAMA 5-1 Cálculo de aproximaciones por diferencias

#### A) Explicaciones

El programa encuentra la aproximación por diferencias para la derivada del orden deseado utilizando los puntos de la retícula especificados por el usuario.

El programa pregunta por: 1) el número de puntos en la retícula que se usarán en la fórmula para la aproximación por diferencias (se puede utilizar un máximo de 10 puntos); 2) los índices de los puntos de la retícula, y 3) el orden de la derivada que se aproxima. El programa supondrá que el intervalo de separación es $h$, sin que el usuario especifique su valor numérico.

Para especificar la aproximación por diferencias deseada, denotamos el número de puntos en la retícula por $L$, el orden de la derivada por $p$ y las abscisas de los puntos de la retícula como $x_\alpha, x_\beta, \dots, x_\lambda$, donde $x_\alpha = \alpha h, x_\beta = \beta h, \dots, x_\lambda = \lambda h$.

El algoritmo funciona incluso cuando $\alpha, \beta, \dots, \lambda$ no son enteros sino que representan cualquier valor decimal positivo o negativo. Las coordenadas de los puntos en la retícula se hacen $x_i = ih$, $i = \alpha, \beta, \dots, \lambda$ y la derivada del orden especificado se evalúa en $x = 0$. Por ejemplo, si se va a evaluar $f''(0)$ utilizando puntos en $x = -2$, $0.5$ y $1.5$, entonces hacemos $L = 3$, $\alpha = -2$, $\beta = 0.5$ y $\gamma = 1.5$ con $h = 1$.

La salida del programa está dada en la forma de los coeficientes de la ecuación (5.3.1); es decir, $a_\alpha, a_\beta, \dots, a_\lambda$, y $c_1$ y $c_2$ de la ecuación (5.3.2). Como se vio en la sección 5.3, el segundo término de la ecuación (5.3.2) se debe ignorar si $c_1 = 0$.

Los coeficientes de las ecuaciones lineales se guardan en el arreglo `A(K, L)` al igual que en `B(K, L)`. El primero se utiliza para la solución, mientras que el segundo se reserva para su uso posterior. Las ecuaciones lineales se resuelven mediante la subrutina de la eliminación gaussiana. En el capítulo 6 se explicarán más detalles del esquema de eliminación gaussiana. Al regresar de la subrutina, la solución de las ecuaciones lineales se guarda en `A(K, KM + 1)`, K = 1, 2, ..., KM. Los coeficientes de la aproximación por diferencias son números decimales en primera instancia. Para expresarlos en forma racional con coeficientes enteros en el numerador y denominador, se hacen algunos cálculos adicionales.

#### B) Variables

* **KM:** número de puntos en la retícula que se utilizará en la aproximación por diferencias ($L$)
* **EL(K):** índice de los puntos de la retícula para el K-ésimo punto contando desde la izquierda (valores de $\alpha, \beta, \gamma, \dots$)
* **DR:** orden de la derivada que se aproximará ($p$)
* **A(K, L):** coeficientes de la ecuación lineal [véase la ecuación (5.3.5)]
* **C(K):** coeficientes del K-ésimo valor de la función en el numerador de la aproximación por diferencias para $K \le KM$: coeficiente del término del error para $K > KM$
* **F:** recíproco del denominador de la aproximación por diferencias

#### C) Listado

```fortran
C-----CSL/F5-1    CALCULO DE APROXIMACIONES POR DIFERENCIAS
      COMMON N, A(10,11), EL(10), B(10,11),C(10) , CF(11)
      MT=6
      PRINT *, 'CSL/F5-1    CALCULO DE APROXIMACIONES POR DIFERENCIAS'
      PRINT *
90    PRINT *, '¿NUMERO DE PUNTOS? '
      READ*, KM
92    IF (KM.GT.2.OR.KM.LE.10) GOTO 95
      PRINT *, ' ENTRADA NO VALIDA. POR FAVOR REPITA LA ENTRADA '
      GO TO 90
95    PRINT *, ' ¿NUMERO DE PUNTOS EN LA RETICULA? ',
     1 '  (OPRIMA ENTER DESPUES DE CADA NUMERO)'
      DO  K=1 , KM
        print 98, k
        READ *,  EL(K)
98      FORMAT( ' ¿INDICE DEL PUNTO?, ',I2,'?')
      END DO
103   PRINT *, ' DE EL ORDEN DE LA DERIVADA '
105   READ *, KDR
106   Z=1.0
      DO I =1,KDR
        Z=Z*FLOAT(I)
      END DO
110   DO 130  K=1, KM+2
        DO L=1, KM
          IF(K.EQ.1) A(K,L)=1.0
          IF(K.GT.1) A(K,L)=EL(L)**(K-1)       ! Preparación de los coeficientes de la matriz
          B(K,L)=A(K,L)                        ! Almacenamiento de los mismos en B (1, J)
        END DO
130   CONTINUE
135   FF=1
        DO K=1, KM
          A(K,KM+1)=0
          IF (K-1.EQ.KDR)  A(K,KM+1)=Z         ! Término no homogéneo
        END DO
140   N=KM
160   PRINT *
      KMP2=KM+2
      CALL GAUSS
170   DO 190 K=1, KM+2
        C(K)=0.0
        DO L=1, KM
          C(K)=C(K)+B(K,L)*A(L,KM+1)           ! Coeficientes del término del error
        END DO
190   CONTINUE
191   F=1000.0
      DO 194  K=1, KM
        IF( A(K,KM+1).EQ.0)  GOTO 194
        IF( ABS(A(K,KM+1)).LT.0.0001)    GOTO 194
        U=ABS(A(K,KM+1))
        IF (U .LT.F) F=U
194   CONTINUE
C                                              !-- Coeficientes de la fórmula de diferencias
      DO K=1, KM
        CF(K)=A(K,KM+1)/F
      END DO
198   print 197
197   FORMAT(' ESQUEMA DE DIFERENCIAS ')
      DO 210 K=1,KM
        FINV=1.0/F
        print 7002, CF(K),FINV,KDR,EL(K)
7002  FORMAT(1X,'+[', F10.5,'/(',F8.5,' H**',I1,')] F(',F6.3,'H)')
210   CONTINUE
      print *
      print 7005
7005  FORMAT(' TERMINO DEL ERROR ')
217   DO K=1, KM+2
        IF (ABS(C(K)).LT.0.00000001) C(K)=0
      END DO
      DD=1.0
      DO K=1, KM
        DD=DD*FLOAT(K)
      END DO
C
      DO K=KM+1, KM+2
        CM=-C(K)
        CPDD=-C(K)/DD
        KM1=K-1
        NH=KM1-KDR
        IF(K.EQ.KM+1.AND.CM.NE.0)
     &            print 7020, CM,DD, NH, KM1     ! Imprime los términos del error
        IF(K.EQ.KM+2) print 7020, CM,DD, NH, KM1
        DD=DD*FLOAT(K)
      END DO
7020  FORMAT( 5X,'(', F10.5,'/',F10.5,')H**',I1,2X,'F^(',  I1, ')')
      print 7030
7030  FORMAT( /'--------------------------------------------------')
      PRINT *
      PRINT *,' OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR '
      READ *,KK
      IF (KK.EQ.1) GO TO 90
      PRINT *
      END
C***********************
      SUBROUTINE GAUSS        !-- Solución de ecuaciones simultáneas
      COMMON N, A(10,11), EL(10), B(10,11),C(10)
      NM=N-1
      N1=N+1
      DO 1085 I=1,NM
        IPV=I
        I1=I+1
        DO J=I1, N
          IF(ABS(A(IPV,I)).LT.ABS(A(J,I))) IPV=J
        END DO
        IF (IPV.EQ.I) GOTO 1060
        DO JC=1,N1
          TM=A(I,JC)
          A(I,JC)=A(IPV,JC)
          A(IPV,JC)=TM
        END DO
1060    DO 1080 JR=I1, N                       ! Comienza la eliminación hacia adelante
          IF (A(JR,I).EQ.0) GOTO 1080
          IF(A(I,I).EQ.0) PRINT *,I,JR,(A(I,JJJ),JJJ=1,3)
          R=A(JR,I)/A(I,I)
          N1=N+1
1075      DO KC=I1, N1
            A(JR,KC)=A(JR,KC)-R*A(I,KC)
          END DO
1080    CONTINUE
1085  CONTINUE
1090  IF (A(N,N).EQ.0) GOTO 1200
1100  A(N,N+1)=A(N,N+1)/A(N,N)                 ! Comienza la sustitución hacia atrás
1110  DO 1130 NV=NM, 1, -1
        VA=A(NV,N+1)
        NV1=NV+1
        DO K=NV1, N
          VA=VA-A(NV,K)*A(K,N+1)
        END DO
        A(NV,N+1)=VA/A(NV,NV)
1130  CONTINUE
      RETURN
1200  print 1210
1210  FORMAT( ' LA MATRIZ ES SINGULAR')
      STOP
      END

```

#### D) Ejemplo de salida

```text
CSL/F5-1    CALCULO DE APROXIMACIONES POR DIFERENCIAS

¿NUMERO DE PUNTOS?
3
¿NUMERO DE PUNTOS EN LA RETICULA? (OPRIMA ENTER DESPUES DE CADA NUMERO)
¿INDICE DEL PUNTO? 1?
0
¿INDICE DEL PUNTO, 2?
1
¿INDICE DEL PUNTO, 3?
2
 DE EL ORDEN DE LA DERIVADA
1

ESQUEMA DE DIFERENCIAS
+[  -3.00000/( 2.00000 H**1)] F( 0.000H)
+[   4.00000/( 2.00000 H**1)] F( 1.000H)
+[  -1.00000/( 2.00000 H**1)] F( 2.000H)

TERMINO DE ERROR
     (   2.00000/    6.00000)H**2  F^(3)
     (   6.00000/   24.00000)H**3  F^(4)
--------------------------------------------------
OPRIMA 1 PARA CONTINUAR, O 0 PARA TERMINAR
¿NUMERO DE PUNTOS?
5
¿NUMERO DE PUNTOS EN LA RETICULA?  (OPRIMA ENTER DESPUES DE CADA NUMERO)
¿INDICE DEL PUNTO?, 1?
-2
¿INDICE DEL PUNTO?, 2?
-1
¿INDICE DEL PUNTO?, 3?
0
¿INDICE DEL PUNTO?, 4?
1
¿INDICE DEL PUNTO?, 5?
2
 DE EL ORDEN DE LA DERIVADA
2

ESQUEMA DE DIFERENCIAS
+[  -1.00000/(12.00000 H**2)] F(-2.000H)
+[  16.00000/(12.00000 H**2)] F(-1.000H)
+[ -30.00000/(12.00000 H**2)] F( 0.000H)
+[  16.00000/(12.00000 H**2)] F( 1.000H)
+[  -1.00000/(12.00000 H**2)] F( 2.000H)

TERMINO DE ERROR
     (   0.00000/  120.00000)H**3  F^(5)
     (   8.00000/  720.00000)H**4  F^(6)
--------------------------------------------------
OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR
