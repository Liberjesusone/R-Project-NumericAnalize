# PROGRAMA 7-1 Método de interpolación

## A) Explicaciones

Este programa transforma la función característica de una matriz en una serie de potencias mediante el método de interpolación.

Antes de correr el programa, hay que definir en los enunciados `DATA` el orden de la matriz, $N$, así como los elementos de la matriz. Cuando se ejecuta el programa, se le da a la computadora en forma interactiva, el tamaño $\lambda$ del intervalo. Después de esto, se calculan los $N + 1$ valores del determinante mediante la subrutina `DETMNT` y se almacenan en el arreglo $\text{FF}(JJ)$. A continuación se calcula la tabla de diferencias y se almacena en $\text{DF}(I, J)$. Los coeficientes del polinomio se calculan mediante los coeficientes de Markov, los que a su vez se calculan internamente cuando se necesitan. La precisión de los coeficientes del polinomio puede verse afectada por la entrada para el intervalo de $\lambda$.

## B) Variables

* $A(I, J)$: elementos de la matriz
* $N$: orden de la matriz
* $\text{DE}$: incremento de $\lambda$
* $\text{RA}(I)$: valores discretos de $\lambda$ para los que se calcula $f(\lambda)$
* $\text{DF}(I, J)$: tabla de diferencias
* $\text{CC}(I)$: coeficiente de la potencia
* $\text{MV}(L)$: coeficiente de Markov

## C) Listado

```fortran
C-----CSL/F7-1.FOR      METODO DE INTERPOLACION (FORTRAN)
C
      DIMENSION A(20,20),FF(0:20),RA(0:20),DF(0:20,0:20),B(20,20)
      DIMENSION CC(0:20)
      REAL MV(0:20)
      PRINT *
      PRINT *, 'CSL/F7-1       METODO DE INTERPOLACION (FORTRAN)'
      PRINT *
C-----DEFINICION DE LA MATRIZ
      DATA N/4/                   ! Defina el orden de la matriz
      DATA (B(1,J),J=1,4)/4,3,2,1/
      DATA (B(2,J),J=1,4)/3,3,2,1/
      DATA (B(3,J),J=1,4)/2,2,2,1/
      DATA (B(4,J),J=1,4)/1,1,1,1/
      PRINT *, 'ORDEN DE LA MATRIZ =',N
      PRINT *, 'MATRIZ '
      DO I=1,N
        PRINT 56, (B(I,J),J=1,N)
      END DO
      PRINT *
      PRINT *, 'DELTA LAMBDA ?'
      READ *,DE                   ! Incremento de lambda para la tabla de diferencias
      DO JJ=0,N
        RA(JJ)=JJ*DE
        PRINT *
        PRINT *, 'LAMBDA= ',RA(JJ)
        DO J=1,N
          DO I=1,N
            A(I,J)=B(I,J)
          END DO
          A(J,J)=A(J,J)-RA(JJ)
        END DO
        CALL DETMNT(N,A,S)
        FF(JJ)=S
      END DO
   56 FORMAT(1X, 1P6E12.5)
C--La siguiente parte calcula la tabla de diferencias hacia adelante
      DO I=0,N
        DF(I,0)=FF(I)             ! Inicialización de la tabla de diferencias
      END DO
      M=N
      DO J=1,N
        M=M-1
        DO I=0,M
          DF(I,J)=DF(I+1,J-1)-DF(I,J-1) ! Tabla de diferencias
        END DO
      END DO
      PRINT *
      PRINT *, ' TABLA DE DIFERENCIAS PARA LOS DETERMINANTES '
      DO I=0,N
        WRITE (*,'(F8.4, 1P6E11.3)') RA(I), (DF(I,J),J=0,N-I)
      END DO
C--La siguiente rutina es para calcular los coeficientes de las potencias por medio de los coeficientes de Markov
      PRINT *
      PRINT *, ' COEFICIENTES DE MARKOV '
      DO I=0,N
        CC(I)=0
        MV(I)=0
      END DO
      MV(1)=1                     ! Se inicializan los coeficientes de Markov
      CC(0)=DF(0,0)
      CC(1)=DF(0,1)
      DO K=2,N
        DO L=K,1,-1
          MV(L)=(MV(L-1)-(K-1)*MV(L))/K
          CC(L)=CC(L)+MV(L)*DF(0,K) ! Coeficientes de Markov
        END DO
      END DO
      PRINT 223, (MV(L),L=1,K)
  223 FORMAT(1X, 6F12.7)
      END DO
      PRINT *
      PRINT *, '--- RESULTADO FINAL ---'
      PRINT *, 'POTENCIA,    COEFICIENTES'
      DO I=0,N
        CC(I)=CC(I)/DE**I
        WRITE (6,'(2X,I3,5X,F10.4)') I, CC(I)
      END DO
      PRINT *
      PRINT *,'---------------------------------------------------------'
```

### Subrutina DETMNT

```fortran
      SUBROUTINE DETMNT (N,A,DET)   ! Cálculo del determinante
      DIMENSION A(20,20)
      INTEGER PV,PC
      PC=0                          ! Inicialización del contador de pivoteo
      DO I=1,N-1                    ! Comienza la eliminación hacia atrás
        PV=I
        DO J=I+1, N
          IF (ABS(A(PV,I)) .LT. ABS(A(J,I))) PV=J
        END DO
        IF (PV.NE.I) THEN
          DO JC=1,N
            TM=A(I,JC)
            A(I,JC)=A(PV,JC)
            A(PV,JC)=TM
          END DO
          PC=PC+1
        END IF
        DO JR=I+1, N
          IF (A(JR,I) .NE. 0) THEN
            R=A(JR,I)/A(I,I)
            DO KC=I+1,N
              A(JR,KC)=A(JR,KC)-R*A(I,KC)
            END DO
          END IF
        END DO
      END DO                        ! Fin de la eliminación hacia atrás
      IF (A(N,N) .EQ. 0) GOTO 1200
      DET=1                         ! Inicialización del determinante
      DO I=1,N
        DET=DET*A(I,I)
      END DO
      IF (PC .NE. INT(PC/2)*2) DET=-DET
      PRINT *, ' DETERMINANTE = ',DET
      PRINT *, ' NUMERO DE PIVOTEOS= ',PC
      RETURN
 1200 PRINT *, ' LA MATRIZ ES SINGULAR'
      RETURN
      END
```

## D) Ejemplo de salida

```text
CSL/F7-1       METODO DE INTERPOLACION (FORTRAN)

ORDEN DE LA MATRIZ =        4
MATRIZ 
 4.00000E+00  3.00000E+00  2.00000E+00  1.00000E+00
 3.00000E+00  3.00000E+00  2.00000E+00  1.00000E+00
 2.00000E+00  2.00000E+00  2.00000E+00  1.00000E+00
 1.00000E+00  1.00000E+00  1.00000E+00  1.00000E+00

DELTA LAMBDA ?
0.5

LAMBDA=    0.000000E+00
DETERMINANTE =        1.000000
NUMERO DE PIVOTEOS =        0

LAMBDA=    0.5000000
DETERMINANTE =        6.2499966E-02
NUMERO DE PIVOTEOS =        1

LAMBDA=    1.000000
DETERMINANTE =       -8.9406981E-08
NUMERO DE PIVOTEOS =        0

LAMBDA=    1.500000
DETERMINANTE =       -4.437500
NUMERO DE PIVOTEOS =        1

LAMBDA=    2.000000
DETERMINANTE =      -17.00000
NUMERO DE PIVOTEOS =        1

TABLA DE DIFERENCIAS PARA LOS DETERMINANTES
 0.0000  1.0000E+00 -9.375E-01  8.750E-01 -5.250E+00  1.500E+00
 0.5000  6.250E-02 -6.250E-02 -4.375E+00 -3.750E+00
 1.0000 -8.941E-08 -4.438E+00 -8.125E+00
 1.5000 -4.438E+00 -1.256E+01
 2.0000 -1.700E+01

COEFICIENTES DE MARKOV
-0.5000000   0.5000000
 0.3333333  -0.5000000   0.1666667
-0.2500000   0.4583333  -0.2500000   0.0416667

--- RESULTADO FINAL ---
POTENCIA,    COEFICIENTES
  0         1.0000
  1        -7.0000
  2        15.0000
  3       -10.0000
  4         1.0000
---------------------------------------------------------
```

---

# PROGRAMA 7-2 Householder/bisección

## A) Explicaciones

Este programa calcula los valores propios de una matriz real simétrica transformando la matriz en una matriz tridiagonal simétrica y calculando en seguida los valores propios de ésta mediante el método de bisección.

El orden y los elementos de la matriz se especifican en los enunciados `DATA`. La reducción de la matriz propuesta a la forma tridiagonal se lleva a cabo en el programa principal por el esquema de Householder. Los enunciados en este proceso se pueden comparar fácilmente con las ecuaciones en el texto. Aquí, la matriz tridiagonal es un caso especial de la matriz de Hessenberg. Después de esto, se calculan los valores propios de la matriz tridiagonal. El determinante de una matriz tridiagonal se calcula utilizando la ecuación (7.3.15) en la subrutina `BISEC`. Los resultados finales se imprimen.

## B) Variables

* $N$: orden de la matriz
* $A(I, J)$: elementos de la matriz
* $S$: $s$
* $\text{SSR}$: $\sqrt{s}$
* $\text{UAU}$: $s = u^T A u$
* $T(I, J)$: memoria temporal de matrices
* $\text{EI}(L, J)$: J-ésimo valor propio de la L-ésima matriz tridiagonal principal
* $XL$: cota inferior para las estimaciones del valor propio
* $XH$: cota superior para las estimaciones del valor propio

## C) Listado

```fortran
C-----CSL/F7-2.FOR      HOUSEHOLDER/TRIDIAGONAL '
      DIMENSION A(0:10,0:10),U(0:10),T(0:10,0:10)
      DIMENSION G(0:10),EI(0:10,0:10)
      COMMON L
      PRINT *
      PRINT *, 'CSL/F7-2                HOUSEHOLDER/TRIDIAGONAL '
      PRINT *, 'SOLO PARA MATRICES SIMETRICAS ='
      PRINT *
C-----------------------------------------------------------------------
C--ESTE PROGRAMA CALCULA LOS VALORES PROPIOS DE UNA MATRIZ SIMETRICA EN DOS
C--ETAPAS; ETAPA 1:
C--                 REDUCCION DE UNA MATRIZ A UNA FORMA TRIDIAGONAL MEDIANTE EL
C--                 ESQUEMA DE HOUSEHOLDER
C--ETAPA 2:
C--                 CALCULO DE LOS VALORES PROPIOS DE LA MATRIZ TRIDIAGONAL MEDIANTE EL
C--                 ESQUEMA DE BISECCION COMBINADO CON LA INTERPOLACION LINEAL
C-----------------------------------------------------------------------
   20 N=4                           ! Define el orden de la matriz
      DATA (A(1,J),J=1,4)/4, 3, 2, 1/ ! Matriz de N x N
      DATA (A(2,J),J=1,4)/3, 3, 2, 1/
      DATA (A(3,J),J=1,4)/2, 2, 2, 1/
      DATA (A(4,J),J=1,4)/1, 1, 1, 1/
      PRINT *, '----------------  REDUCCION DE UNA MATRIZ A LA FORMA DE HESSENBERG '
      PRINT *, '                 O TRIDIAGONAL MEDIANTE EL ESQUEMA DE HOUSEHOLDER '
      PRINT *
      PRINT *, 'MATRIZ ORIGINAL (DEBE SER SIMETRICA) '
      PRINT *
      DO I=1,N
        PRINT 52, (A(I,J),J=1,N)
      END DO
   52 FORMAT(1X, 1P6E12.5)
      DO IR=1,N-2                   ! Comienza el esquema de Householder.
        S=0
        DO I=1,N
          U(I)=0
          IF (I .GT. IR+1) U(I)=A(I,IR)
          IF (I .GT. IR)   S=S+A(I,IR)*A(I,IR)
        END DO
        W=1
        IF (A(IR+1,IR) .LT. 0) W=-1
        SSR=SQRT(S)
        PRINT *
        H=S+ABS(A(IR+1,IR))*SSR
        U(IR+1)=A(IR+1,IR)+SSR*W
        UAU=0
        DO I=1,N
          DO J=1,N
            UAU=UAU+U(I)*A(I,J)*U(J)
            IF ((I .LE. IR).AND.(J .LE. IR)) THEN
              T(I,J)=A(I,J)
              GOTO 710
            ENDIF
            IF ((J.EQ.IR).AND.(I .GE. IR+2)) THEN
              T(I,J)=0
              GOTO 710
            ENDIF
            B23=0
            DO K=1,N
              B23=B23-(U(I)*A(K,J)+A(I,K)*U(J))*U(K)
            END DO
            T(I,J)=A(I,J)+B23/H
  710     END DO
        END DO
        UAU=UAU/H/H
        DO I=1,N
          DO J=1,N
            A(I,J)=T(I,J)+UAU*U(I)*U(J)
            IF (ABS(A(I,J)) .LT. .000001) A(I,J)=0
          END DO
        END DO
      END DO
      PRINT *, ' MATRIZ DE HESSENBERG O TRIDIAGONAL '
      PRINT *
      DO I=1,N
        PRINT 52, (A(I,J),J=1,4)    ! Imprime la matriz de Hessenberg
      END DO
      PRINT *
      PRINT *, 'PARA CONTINUAR, OPRIMA 1 Y LA TECLA ENTER '
      READ*, DUMM
      PRINT *
      KM=N
      DO L=1,KM
        IF (L.EQ.1) THEN
          EI(1,1)=A(1,1)
        ELSE
          DO J=1,L
            XL=EI(L-1,J-1)
            XH=EI(L-1,J)
            KM=J
            CALL BISEC(G,A,XL,XH,XM)
            EI(L,J)=XM
          END DO
        END IF
        EI(L,0)=-99
        EI(L,L+1)=99
        IF (L.NE.N) THEN
          PRINT *
          PRINT 2130, L,L
 2130     FORMAT(' VALORES PROPIOS DE LA ',I2,' X',I2,
     1           '    SUBRIATRIZ PRINCIPAL   ')
        ELSE
          PRINT *
          PRINT *, ' RESULTADOS FINALES (VALORES PROPIOS DE LA MATRIZ COMPLETA) '
          PRINT *, '---------------------------------------------------------'
        END IF
        PRINT 52, (EI(L,I),I=1,L)
      END DO
      PRINT *, '---------------------------------------------------------'
      PRINT *
      PRINT *
      PRINT *
      STOP
      END
```

### Subrutinas BISEC y DETERM

```fortran
C*******************************************************
      SUBROUTINE BISEC(G,A,XL,XH,XM)     ! Bisección
      DIMENSION A(0:10,0:10), G(0:10)  ! Calcula las raíces de un determinante
      COMMON L
      KA=0
   30 CALL DETERM(G,A,XL,YL)
      CALL DETERM(G,A,XH,YH)
      KA=KA+1
      IF (KA .GT. 99) RETURN
      DX=XH-XL
      IF (DX .LT. .0000001) RETURN
      IF (DX .GT. 1) THEN
        XM=(XL+XH)/2                ! Esquema de bisección
        CALL DETERM(G,A,XM,YM)
        GOTO 30
      END IF
      XB=XM
      XM=(XL*YH-XH*YL)/(YH-YL)        ! Esquema de interpolación lineal
      CALL DETERM(G,A,XM,YM)
      IF (ABS(XB-XM) .LT. .000001) RETURN
      IF (YL*YM .LT. 0) THEN
        XH=XM
        YH=YM
        GOTO 80
      ENDIF
      XL=XM
      YL=YM
      GOTO 80
      END
C*******************************************************
      SUBROUTINE DETERM(G,A,X,SL)       ! Calcula el determinante de una matriz tridiagonal
      DIMENSION A(0:10,0:10),G(0:10)
      COMMON L
      G(0)=1
      IF (L.EQ.1) RETURN
      G(1)=A(1,1)-X
      IF (L.EQ.1) RETURN
      DO K=2,L
        G(K)=(A(K,K)-X)*G(K-1)-A(K,K-1)*A(K,K-1)*G(K-2)
      END DO
      SL=G(L)
      RETURN
      END
```

## D) Ejemplo de salida

```text
CSL/F7-2        HOUSEHOLDER/TRIDIAGONAL
(SOLO PARA MATRICES SIMETRICAS)
REDUCCION DE UNA MATRIZ A LA FORMA DE HESSENBERG
O TRIDIAGONAL MEDIANTE EL ESQUEMA DE HOUSEHOLDER

MATRIZ ORIGINAL (DEBE SER SIMETRICA)
 4.00000E+00  3.00000E+00  2.00000E+00  1.00000E+00
 3.00000E+00  3.00000E+00  2.00000E+00  1.00000E+00
 2.00000E+00  2.00000E+00  2.00000E+00  1.00000E+00
 1.00000E+00  1.00000E+00  1.00000E+00  1.00000E+00

MATRIZ DE HESSENBERG O TRIDIAGONAL
 4.00000E+00 -3.74166E+00  0.00000E+00  0.00000E+00
-3.74166E+00  5.00000E+00  4.62911E-01  0.00000E+00
 0.00000E+00  4.62911E-01  6.66666E-01 -8.90870E-02
 0.00000E+00  0.00000E+00 -8.90869E-02  3.33333E-01

VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 1 x 1
 4.00000E+00

VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 2 x 2
 7.25082E-01  8.27492E+00

VALORES PROPIOS DE LA SUBMATRIZ PRINCIPAL DE 3 x 3
 3.81085E-01  9.94722E-01  8.29086E+00

RESULTADOS FINALES (VALORES PROPIOS DE LA MATRIZ COMPLETA)
---------------------------------------------------------
 2.83250E-01  4.26021E-01  9.99996E-01  8.29086E+00
---------------------------------------------------------