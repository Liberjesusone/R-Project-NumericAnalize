# PROGRAMAS

## PROGRAMA 4-1 Reglas extendidas del trapecio y de Simpson

### A) Explicaciones

El PROGRAMA 4-1 integra una función analítica, ya sea mediante la regla extendida del trapecio o mediante la regla extendida de Simpson, según la elección del usuario. Antes de correr el programa, el usuario debe definir el integrando en el subprograma FUNC. El usuario puede dar como entrada la elección de un método de integración, los límites de integración y el número de intervalos, en forma interactiva desde el teclado. Si el número de intervalos para la regla de Simpson es impar, se utiliza la regla de 3/8 para los primeros tres intervalos de la retícula y después se utiliza la regla extendida de 1/3 para el resto del dominio.

Cuando se corre el programa, la computadora envía un mensaje para recordar al usuario que la función a integrar debe definirse en el subprograma FUNC, después le pregunta el método que usará para integrar; si es el de SIMPSON la entrada debe ser 1. A continuación, el programa pregunta los valores de N, A y B.

### B) Variables

* **ISIMP:** especificación del método
* ISIMP = 0 (Regla del trapecio)
* ISIMP = 1 (Regla de Simpson)


* **A, B:** límite inferior y superior de integración, respectivamente
* **N:** número de intervalos en la retícula
* **H:** espaciamiento, H = (B - A)/N
* **W:** valores de los pesos en las fórmulas de integración
* **S, SS:** integral
* **II:** último punto de la retícula para la regla de 3/8

### C) Listado

```fortran
C----CSL/F4-1.FOR   REGLAS DEL TRAPECIO Y DE SIMPSON
      COMMON A,B,H
      CHARACTER SIMP*6
      PRINT *
      PRINT *,'CSL/F4-1   REGLAS DEL TRAPECIO Y DE SIMPSON'
      PRINT *
      PRINT *,'LA FUNCION A INTEGRAR SE DEBE CODIFICAR EN
      PRINT *,' EL SUBPROGRAMA LLAMADO FUNC
      PRINT *
10    PRINT *,'OPRIMA 0 PARA EL TRAPECIO, 1 PARA SIMPSON
      READ *, ISIMP
      PRINT *,'¿NUMERO DE INTERVALOS?
      READ *,N
135   IF (N .GT.0 .AND. ISIMP.EQ.0) GOTO 140
      IF (ISIMP.EQ.1 .AND. N.GT.1) GOTO 140
      PRINT *,'LA ENTRADA NO ES VALIDA, REPITA '
      GO TO 10
140   PRINT *,'¿LIMITE INFERIOR DE INTEGRACION?      '
      READ *,A
150   PRINT *,'LIMITE SUPERIOR DE INTEGRACION ? '
      READ *,B
160   H=(B-A)/N
      IF (ISIMP.EQ.0) THEN
         CALL  TRAPZ(S,N)       !-- Se eligió la regla del trapecio.
         GOTO 200
      ELSE
         CALL  SIMPS(S,N)       !-- Se eligió la regla de Simpson.
      END IF
200   PRINT *,'--------------------------------------------------'
210   PRINT *,'RESULTADO FINAL=',S
220   PRINT *,'--------------------------------------------------'
      PRINT *
      PRINT*
      PRINT*,'OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR '
      READ *,K
      IF(K.EQ.1) GOTO 10
      PRINT*
      END
C**********************************
      SUBROUTINE TRAPZ(S,N)     !-- Regla del trapecio
      COMMON A,B,H
      S=0
      DO 10 I=0,N
        X=A+I*H
        W=2
        IF(I.EQ.0 .OR. I.EQ.N) W = 1
        S=S+W*FUNC(X)
c       PRINT *,I,X,H,FUNC(X),W
10    CONTINUE
      S=S*H/2
      RETURN
      END
C**********************************
      SUBROUTINE SIMPS(SS,N)
      COMMON A,B,H
      S=0
      SS=0
      IF (N/2*2.EQ.N) THEN
         LS = 0
         GOTO 35
      END IF
      LS=3
      DO 30 I=0,3             !  Regla de 3/8 de Simpson si N es impar
         X=A+H*I
         W=3
         IF (I.EQ.0 .OR. I.EQ.3) W=1
         SS=SS+W*FUNC(X)
30    CONTINUE
      SS=SS*H*3/8
      IF (N.EQ.3) RETURN
35    DO 40 I=0, N-LS         !  Regla de 1/3 de Simpson
         X=A+H*(I+LS)
         W=2
         IF (INT(I/2)*2+1.EQ.I) W=4
         IF (I.EQ.0 .OR. I.EQ.N-LS)  W=1
         S=S+W*FUNC(X)
40    CONTINUE
      SS=SS+S*H/3
      RETURN
      END
C**********************************
      FUNCTION FUNC(X)        !-- Evalúa la función a integrar
      FUNC = (1 + (X/2)**2)**2*3.14159
      RETURN
      END

```

### D) Ejemplo de salida

```text
CSL/F4-1       REGLAS DEL TRAPECIO Y DE SIMPSON
LA FUNCION A INTEGRAR SE DEBE CODIFICAR EN EL SUBPROGRAMA LLAMADO FUNC
OPRIMA 0 PARA EL TRAPECIO, 1 PARA SIMPSON
0
¿NUMERO DE INTERVALOS?
10
¿LIMITE INFERIOR DE INTEGRACION?
0
¿LIMITE SUPERIOR DE INTEGRACION?
2
--------------------------------------------------
RESULTADO FINAL    11.77047
--------------------------------------------------
OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR
1
OPRIMA 0 PARA EL TRAPECIO, 1 PARA SIMPSON
1
¿NUMERO DE DATOS?
5
LIMITE INFERIOR DE INTEGRACION ?
0
LIMITE SUPERIOR DE INTEGRACION ?
2
--------------------------------------------------
RESULTADO FINAL    11.73095
--------------------------------------------------
OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR

```

Aquí tienes la transcripción de las nuevas imágenes al formato Markdown (.md), manteniendo la misma estructura y formato de los bloques de código:

---

## PROGRAMA 4-2 Fórmulas cerradas de Newton-Cotes

### A) Explicaciones

El PROGRAMA 4-2 lleva a cabo la integración numérica utilizando las fórmulas cerradas de Newton-Cotes, mientras que el PROGRAMA 4-3 utiliza las fórmulas abiertas.
Antes de correr el programa, el usuario debe definir la función a integrar en el subprograma FUN, donde aparece F = SIN (X) como ejemplo. Las instrucciones DATA contienen los valores de W, Q y R copiados de la tabla 4.2. El programa pregunta por el orden del método y los límites de integración.

### B) Variables

* **N:** orden de la fórmula de Newton-Cotes, $2 \le N \le 10$
* **A y B:** límites inferior y superior de la integral
* **W(J, N):** J-ésimo factor de peso en la fórmula de orden N (tablas 4.2 y 4.3)
* **Q, R:** numerador y denominador de $\alpha = Q/R$ (tablas 4.2 y 4.4)
* **H:** espaciamiento, (B - A)/N
* **I:** integral (respuesta final)

### C) Listado

```fortran
C----CSL/F4-2.FOR        FORMULA CERRADA DE NEWTON-COTES
      DIMENSION W(0:20,10)
      PRINT *
      PRINT *,'CSL/F4-2      FORMULA CERRADA DE NEWTON COTES '
      PRINT *
      DATA (W(I,1),I=0,3)/1,1,1,2/
      DATA (W(I,2),I=0,4)/1,4,1,1,3/
      DATA (W(I,3),I=0,5)/1,3,3,1,3,8/
      DATA (W(I,4),I=0,6)/7,32,12,32,7,2,45/
      DATA (W(I,5),I=0,7)/19,75,50,50,75,19,5,288/
      DATA (W(I,6),I=0,8)/41,216,27,272,27,216,41,1,140/
      DATA (W(I,7),I=0,7)/751,3577,1323,2989,2989,1323,3577,751/
      DATA (W(I,7),I=8,9)/7,17280/
      DATA (W(I,8),I=0,8)/989,5888,-928, 10946 ,-4540, 10946 ,-928,5888,989/
      DATA (W(I,8),I=9,10)/4,14175/
      DATA (W(I,9),I=0,7)/2857,15741,1080,19344,5788,5788,19344,1080/
      DATA (W(I,9),I=8,11)/15741,2857,9,89600/
      DATA (W(I,10),I=0,5)/16067, 106300,-48525,272400,-260550, 427368 /
      DATA (W(I,10),I=6,10)/-260550,272400,-48525,106300,16067/
      DATA (W(I,10),I=11,12)/5,299376/
1     PRINT *,'¿NUMERO DE DATOS?          (2 - 10) '
      READ *,K
      PRINT *,'A (LIMITE INFERIOR DE INTEGRACION) ? '
      READ *,A
      PRINT *,'B (LIMITE SUPERIOR DE INTEGRACION) ? '
      READ *,B
      PRINT *
      N=K-1
      Q=W(N+1,N)
      R=W(N+2,N)
      PRINT 10,Q,R
10    FORMAT(' Q=',1PE14.6, '    R=', 1PE14.6)
      PRINT *
      AL=Q/R                  ! alfa
      H=(B-A)/N               ! Intervalo de la retícula
      PRINT*,'--------------------------------------------------'
      PRINT*,'      N         X             F(X)            W     '
      PRINT*,'--------------------------------------------------'
220   ANS=0                   ! Inicialización de la fórmula de Newton-Cotes
      DO 240 J=0, N                           ! - J es el índice de los puntos de la retícula
        X=A+J*H
        F=FUN(X)
        PRINT 250, J,X,F,W(J,N)
        ANS=ANS+F*W(J,N)                !-- Fórmula de Newton-Cotes
240   CONTINUE
250   FORMAT(1X,I5,1P2E15.6, F13.4)
      ANS=ANS*H*AL
      PRINT *
      PRINT*,' ---------------------------------------------------------------------- '
      PRINT*,'       RESULTADO FINAL         I=' , ANS
      PRINT*,' ---------------------------------------------------------------------- '
      PRINT*
      PRINT*,' OPRIMA 1 PARA CONTINUAR, O 0 PARA TERMINAR '
      READ *,K
      IF(K.EQ.1) GOTO 1
      PRINT*
      END
C**********************************
      FUNCTION FUN(X)     ! -- Evalúa la función a integrar
      FUN=SIN(X)
      RETURN
      END

```

### D) Ejemplo de salida

```text
CSL/F4-2      FORMULA CERRADA DE NEWTON-COTES

¿ NUMERO DE DATOS ? (2 - 10)
6
A (LIMITE INFERIOR DE INTEGRACION) ?
0
B (LIMITE SUPERIOR DE INTEGRACION) ?
2
Q=   5.000000E+00      R=   2.880000E+02
--------------------------------------------------
      N         X             F(X)            W     
--------------------------------------------------
    0    0.000000E+00    0.000000E+00      19.0000
    1    4.000000E-01    3.894183E-01      75.0000
    2    8.000000E-01    7.173561E-01      50.0000
    3    1.200000E+00    9.320391E-01      50.0000
    4    1.600000E+00    9.995736E-01      75.0000
    5    2.000000E+00    9.092974E-01      19.0000
--------------------------------------------------
      RESULTADO FINAL         I=   1.416117
--------------------------------------------------
OPRIMA 1 PARA CONTINUAR O 0 PARA TERMINAR

```