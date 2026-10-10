# PROGRAMA 6-2 Inversión de una matriz

## A) Explicaciones

Este programa obtiene la inversa de una matriz por medio de la eliminación de Gauss, como se explicó al final de la sección 6.6. Los primeros enunciados `DATA` definen el valor de $N$, que por el momento es 3, con el fin de ejemplificar. El arreglo aumentado también se especifica en un enunciado `DATA`. La inversión de la matriz se lleva a cabo en la subrutina, la cual es esencialmente igual a la eliminación de Gauss del **PROGRAMA 6-1**, excepto que se añade la parte aumentada y que la sustitución hacia atrás se aplica a las $N$ columnas del lado derecho en la matriz aumentada. Cuando se concluye la subrutina, la matriz inversa está incluida en las últimas $N$ columnas del arreglo aumentado $A$. Después de regresar al programa principal, se imprime la matriz inversa.

## B) Variables

* **EPS:** épsilon de la máquina
* **N:** orden de la matriz
* **A(I, J):** arreglo aumentado: las primeras $N$ columnas son la matriz dada y las últimas $N$ columnas son inicialmente la matriz identidad

## C) Listado

```fortran
C-----CSL/F6-2.FOR      INVERSION DE UNA MATRIZ
      DIMENSION A(0:10,0:20)
      PRINT *
      PRINT *, 'CSL/F6-2      INVERSION DE UNA MATRIZ '
      DATA N/3/
      DATA (A(1,J),J=1,6)/2, 1,-3, 1, 0, 0/
      DATA (A(2,J),J=1,6)/-1, 3, 2, 0, 1, 0/
      DATA (A(3,J),J=1,6)/3, 1,-3, 0, 0, 1/
      PRINT *
      PRINT *, 'MATRIZ ORIGINAL '
      PRINT *
      DO I=1, N
        PRINT 20, (A(I,J),J=1,3)
      END DO
 20   FORMAT(1X,1P5E12.5)
      PRINT *
      CALL GAUSS(N,A)
      PRINT *
      PRINT *, ' MATRIZ INVERSA '
      PRINT
      DO I=1,N
        PRINT 20,(A(I,J),J=N+1,N*2)
      END DO
      PRINT *
      STOP
      END
C***********************************
      SUBROUTINE GAUSS(N,A)             ! Subrutina para la eliminación de Gauss
      DIMENSION A(0:10,0:20)
      INTEGER PV
      EPS=1.0
 10   IF(1.0+EPS.GT.1.0) THEN
        EPS=EPS/2.0
        GOTO 10
      END IF
      EPS=EPS*2
      PRINT *, ' EPSILON DE LA MAQUINA =', EPS
      EPS2=EPS*2
      DET=1.0                            ! Inicialización del determinante
      DO 1010 I=1,N-1
        PV=I
        DO J=I+1,N
          IF (ABS(A(PV,I)) .LT. ABS(A(J,I))) PV=J
        END DO
        IF (PV.NE.I) THEN
          DO JC=1,N*2
            TM=A(I,JC)
            A(I,JC)=A(PV,JC)
            A(PV,JC)=TM
          END DO
          DET=-DET
        END IF
        IF (A(I,I).EQ.0) GOTO 1200
C ELIMINACION POR DEBAJO DE LA DIAGONAL ---
        DO JR=I+1,N
          IF (A(JR,I).NE.0) THEN
            R=A(JR,I)/A(I,I)
            DO KC=I+1,N*2
              TEMP=A(JR,KC)
              A(JR,KC)=A(JR,KC)-R*A(I,KC)
              IF(ABS(A(JR,KC)).LT.EPS2*TEMP) A(JR,KC)=0.0
            END DO
          END IF
        END DO
 1010 CONTINUE
      DO I=1,N
        DET=DET*A(I,I)
      END DO
      PRINT *
      PRINT *, ' DETERMINANTE=', DET
      PRINT *
C SUSTITUCION HACIA ATRAS
      IF (A(N,N).EQ.0) GOTO 1200
      DO 1100 M=N+1,N*2
        A(N,M)=A(N,M)/A(N,N)
        DO NV=N-1,1,-1
          VA=A(NV,M)
          DO K=NV+1,N
            VA=VA-A(NV,K)*A(K,M)
          END DO
          A(NV,M)=VA/A(NV,NV)
        END DO
 1100 CONTINUE
      RETURN
 1200 PRINT *, 'LA MATRIZ ES SINGULAR '
      END
```

## D) Ejemplo de salida

```text
CSL/F6-2      INVERSION DE UNA MATRIZ

MATRIZ ORIGINAL
 2.00000E+00 1.00000E+00-3.00000E+00
-1.00000E+00 3.00000E+00 2.00000E+00
 3.00000E+00 1.00000E+00-3.00000E+00
            EPSILON DE LA MAQUINA =   5.9604645E-08
            DETERMINANTE     = 11.00000

MATRIZ INVERSA
-1.00000E+00-9.93411E-09 1.00000E+00
 2.72727E-01 2.72727E-01-9.09091E-02
-9.09091E-01 9.09091E-02 6.36364E-01
```

---

# PROGRAMA 6-3 Descomposición LU

## A) Explicaciones

Este programa descompone una matriz $A$ en la forma $LU$ con pivoteo. En el programa, las matrices $L$ y $U$ se guardan en la forma combinada y se expresan mediante el arreglo $L(I, J)$; vea el párrafo posterior a la ecuación (6.7.12) para una mayor explicación. Sin embargo, para simplificar la demostración, la matriz original $A$ y $L$ se expresan en forma aparte. (Cualquiera de estos arreglos se puede eliminar, como se verá después.)

Para reducir las necesidades de espacio de memoria, se puede eliminar fácilmente uno de los arreglos $A$ o $L$. Por ejemplo, si sólo se va a usar $A(I, J)$, se puede borrar "$L(30, 30)$" de los enunciados de dimensión y remplazar cada uno de los "$L(I, J)$" por "$A(I, J)$" en todo el programa.

Este programa también calcula el determinante de la matriz original con el algoritmo explicado en la sección 6.8.

## B) Variables

* **N:** orden de la matriz
* **A(I, J):** elementos de la matriz $A$
* **EL(I, J):** matriz triangular inferior $L$ y matriz triangular superior $U$ en la forma combinada
* **IP(I):** permutación debida al pivoteo
* **IPC:** número de pivoteos realizados
* **DE:** determinante de la matriz original

## C) Listado

```fortran
C-----CSL/F6-3.FOR      DESCOMPOSICION LU CON PIVOTEO
      DIMENSION A(20,20),EL(20,20),IP(20)
      REAL*4 L
 5    PRINT *
      PRINT *,'CSL/F6-3      DESCOMPOSICION LU '
      PRINT *
      DATA N/3/
      DATA (A(1,J),J=1,3)/ 2, 1, -3/
      DATA (A(2,J),J=1,3)/-1, 3,  2/
      DATA (A(3,J),J=1,3)/ 3, 1, -3/
      PRINT *,' N = ',N
      PRINT *,' MATRIZ ORIGINAL-----------------------------------'
      DO I=1,N
        PRINT 10,(A(I,J),J=1,N)
      END DO
 10   FORMAT(1X, 1P6E12.5)
      PRINT *
      PRINT *,' SI SE DESEA EL PIVOTEO, OPRIMA 1; DE LO CONTRARIO OPRIMA 0 '
      READ *,NP
C--INICIALIZACION DEL PIVOTEO Y LA MATRIZ EL
      IPC=1
      DO I=1,N
        IP(I)=I
        DO J=1,N
          EL(I,J)=0
        END DO
      END DO
      J=1
      IF(NP.EQ.1) CALL PIVOT(N,A,EL,J,IP,IPC)
C  PRIMER RENGLON
      DO J=1,N
        EL(1,J)=A(1,J)
      END DO
C  PRIMERA COLUMNA
      DO I=2,N
        EL(I,1)=A(I,1)/EL(1,1)
      END DO
      DO 80 M=2,N
        IF(NP.EQ.1) CALL PIVOT(N,A,EL,M,IP,IPC)
C                                            ! M-ESIMO RENGLON
        DO J=M,N
          S=0
          DO K=1,M-1
            S=S+EL(M,K)*EL(K,J)
          END DO
          EL(M,J)=A(M,J)-S
        END DO
C                                            ! M-ESIMA COLUMNA
        DO I=M+1,N
          S=0
          DO K=1,M-1
            S=S+EL(I,K)*EL(K,M)
          END DO
          EL(I,M)=(A(I,M)-S)/EL(M,M)
        END DO
 80   CONTINUE
      PRINT *,' PERMUTACION '
      PRINT 361,(IP(I),I=1,N)
361   FORMAT(1X,10I3)
368   PRINT *
      PRINT *,' MATRICES LU EN FORMA COMPACTA '
      DO I=1,N
        PRINT 10,(EL(I,M),M=1,N)
      END DO

      PRINT *
C---CALCULO DEL DETERMINANTE
      DE=1
      DO I=1,N
        DE=DE*EL(I,I)
      END DO
      IF (IPC.EQ.INT(IPC/2)*2) DE=-DE
      PRINT *
      PRINT *,' DETERMINANTE = ',DE
395   PRINT *
      PRINT *,'--------------------------------------------------'
      PRINT *
      PRINT*
      END
C*****************************************
      SUBROUTINE PIVOT(N,A,EL,J,IP,IPC)
      DIMENSION A(20,20),EL(20,20),IP(20)
 405  T=0
      DO 420 K=J,N
        IF (ABS(A(K,J)) .LE.T) GOTO 420
        JJ=K
        T=ABS(A(K,J))
 420  CONTINUE
      IF (JJ.EQ.J) RETURN
 425  IPC=IPC+1
      DO 430 M=1,N
        T=A(J,M)
        A(J,M)=A(JJ,M)
        A(JJ,M)=T
        T=EL(J,M)
        EL(J,M)=EL(JJ,M)
        EL(JJ,M)=T
 430  CONTINUE
      IT=IP(J)
      IP(J)=IP(JJ)
      IP(JJ)=IT
      PRINT *,' NUMERO DE PIVOTEOS =', IPC
      RETURN
      END
```

## D) Ejemplo de salida

```text
CSL/F6-3      DESCOMPOSICION LU

MATRIZ ORIGINAL
 2.00000E+00 1.00000E+00-3.00000E+00
-1.00000E+00 3.00000E+00 2.00000E+00
 3.00000E+00 1.00000E+00-3.00000E+00

SI SE DESEA EL PIVOTEO, OPRIMA 1; DE LO CONTRARIO OPRIMA 0
1
NUMERO DE PIVOTEOS     = 2
PERMUTACION
  3  2  1

MATRICES LU EN FORMA COMPACTA
 3.00000E+00 1.00000E+00-3.00000E+00
-3.33333E-01 3.33333E+00 1.00000E+00
 6.66667E-01 1.00000E-01-1.10000E+00
            DETERMINANTE = 11.00000
```