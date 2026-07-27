# Errores en Métodos Numéricos — Guía de Estudio

> Basado en: Nakamura, S. *Métodos Numéricos Aplicados con Software* y láminas del Prof. Carlos Echeverría, ULA-CeSiMo.

---

## 1. ¿Por qué existen errores en el cálculo numérico?

Los métodos numéricos son **aproximaciones** de procesos matemáticos exactos. Existen tres razones fundamentales por las que siempre habrá error:

1. Las computadoras representan números con **precisión finita** (número fijo de bits).
2. Los algoritmos numéricos aproximan procesos **matemáticamente infinitos** (series, límites, integrales).
3. Los datos de entrada ya contienen **imprecisiones propias** del mundo físico.

---

## 2. Las Fuentes Principales de Error

### 2.1 Error del Modelo
Las hipótesis o condiciones del modelo físico están alejadas de la realidad o son incorrectas. Ejemplo: asumir que una viga es perfectamente rígida cuando en realidad se deforma.

### 2.2 Error del Método
El problema exacto no puede resolverse de forma directa, por lo que se formula una **aproximación**. El error surge de la diferencia entre el modelo exacto y el aproximado.

### 2.3 Error Residual
Originado al truncar **series infinitas** considerando solo una parte finita de ellas. Es la causa del *error de truncamiento*.

### 2.4 Error Inicial
Proviene de los parámetros de entrada cuyos valores son **conocidos de forma aproximada** (datos medidos, constantes físicas redondeadas).

### 2.5 Error de Redondeo
Originado por la **representación finita de los números** en la memoria de la computadora. Todo número real que no sea exactamente representable en binario introduce este error.

### 2.6 Error de Operación
Se divide en dos casos:

- **Operación exacta sobre números aproximados:**
$$\pi - \frac{1}{6} \approx 3.1416 - 0.1667 = 3.2083$$

- **Operación aproximada sobre números exactos:**
$$0.23 \times 10^8 - 1.36 \times 10^{-8} \approx 0.23 \times 10^8$$
El valor verdadero es `22999999.9999999864` — la computadora pierde los decimales.

### 2.7 Error Sistemático
Error que entra de igual modo en cada resultado sin variar las condiciones. Sus causas:
- Defecto del instrumento de medición
- Condiciones del ambiente
- Metodología de la medición
- Precisión limitada del instrumento

### 2.8 Error Causal o Accidental (Fortuito)
Vinculado a factores que sufren pequeñas variaciones **aleatorias** durante el experimento. No es predecible ni controlable.

---

## 3. Las Dos Causas Principales en el Análisis Numérico

El libro de Nakamura y las láminas del profesor enfatizan que en el cálculo numérico computacional, las dos causas dominantes son:

### 3.1 Error de Truncamiento
Se debe a las **aproximaciones utilizadas en la fórmula matemática del modelo**. Ocurre cuando se reemplaza un proceso infinito por uno finito.

**Ejemplo con Series de Taylor:**

La serie exacta es:
$$f(x) = f(a) + h f'(a) + \frac{h^2}{2}f''(a) + \frac{h^3}{6}f'''(a) + \cdots$$

En la práctica se trunca en el término $N$:
$$f(x) \approx f(a) + h f'(a) + \frac{h^2}{2}f''(a) + \cdots + \frac{h^N}{N!}f^{(N)}(a) + O(h^{N+1})$$

Donde el **error de truncamiento** es:
$$O(h^{N+1}) = \frac{h^{N+1}}{(N+1)!} f^{(N+1)}(a + \xi h), \quad 0 \leq \xi \leq 1$$

> **Clave:** Cuanto mayor sea $N$ (más términos), menor es el error de truncamiento, pero mayor es el costo computacional.

### 3.2 Error de Redondeo
Se asocia con el **número limitado de dígitos** con que se representan los números en una computadora.

---

## 4. Representación de Punto Flotante (IEEE 754)

### 4.1 Sistemas de Numeración
Las computadoras trabajan internamente en binario. La equivalencia entre bases:

| Decimal | Binario | Hexadecimal |
|---------|---------|-------------|
| 0–9     | 0–1111  | 0–9, A–F    |

El valor decimal de un número en base $r$: `(abcde.fgh)_r` se calcula como:
$$a \cdot r^4 + b \cdot r^3 + c \cdot r^2 + d \cdot r + e + f \cdot r^{-1} + g \cdot r^{-2} + h \cdot r^{-3}$$

**Ejemplo:** `110101` en binario:
$$1\times2^5 + 1\times2^4 + 0\times2^3 + 1\times2^2 + 0\times2 + 1 = 32+16+4+1 = 53$$

### 4.2 Representación en 64 bits (Doble Precisión)

El estándar IEEE 754 de **doble precisión** usa 64 bits distribuidos así:

```
 63  | 62 ........ 52 | 51 ........................... 0
  S  |   Exponente    |           Mantisa
1 bit|    11 bits     |           52 bits
```

- **Bit de signo (1 bit):** 0 = positivo, 1 = negativo
- **Exponente (11 bits):** sesgado por 1023. Rango: $2^{-1022}$ a $2^{1023}$
- **Mantisa (52 bits):** los dígitos significativos del número

Esto da el rango:
- **Máximo:** $\approx 1.8 \times 10^{308}$
- **Mínimo positivo normalizado:** $\approx 2.2 \times 10^{-308}$
- **Dígitos significativos:** $\approx 15-16$ dígitos decimales

---

## 5. Epsilon de Máquina ($\varepsilon_m$)

**Definición:** El epsilon de máquina es el número más pequeño $\varepsilon$ tal que:
$$1 + \varepsilon > 1 \quad \text{(distinguible de 1 en la computadora)}$$

Es el tamaño del intervalo entre `1` y el siguiente número representable mayor que `1`.

**Código R:**
```r
epsilon <- 1.0
while ((epsilon / 2 + 1) > 1) {
  epsilon <- epsilon / 2
}
cat("Epsilon de maquina calculado:", epsilon, "\n")
cat("Epsilon de maquina (.Machine):", .Machine$double.eps, "\n")
```

**Resultado en doble precisión:**
```
Epsilon de maquina calculado:  2.220446e-16
Epsilon de maquina (.Machine): 2.220446e-16
```

> **Interpretación:** Dos números que difieran en menos de $2.22 \times 10^{-16}$ son **indistinguibles** para la computadora.

---

## 6. Overflow y Underflow

### Overflow (Desbordamiento por arriba)
Ocurre cuando un número **excede el máximo representable**. La computadora devuelve `Inf`.

### Underflow (Desbordamiento por abajo)
Ocurre cuando un número **cae por debajo del mínimo positivo representable**. La computadora lo redondea a `0`.

**Código R:**
```r
# Overflow
x <- 1.0
contador <- 0
while (is.finite(x)) {
  x <- x * 2
  contador <- contador + 1
}
cat("Overflow despues de", contador, "multiplicaciones por 2\n")
cat("Maximo representable:", .Machine$double.xmax, "\n")

# Underflow
y <- 1.0
contador2 <- 0
while (y > 0) {
  y <- y / 2
  contador2 <- contador2 + 1
}
cat("Underflow despues de", contador2, "divisiones por 2\n")
cat("Minimo positivo:", .Machine$double.xmin, "\n")
```

**Resultados de tu máquina:**
```
Overflow despues de 1024 multiplicaciones por 2
Maximo representable: 1.797693e+308

Underflow despues de 1075 divisiones por 2
Minimo positivo: 2.225074e-308
```

> **¿Por qué 1075 y no 1022?** Porque R llega hasta los **números subnormales** antes de llegar a cero. Los subnormales permiten representar números aún más pequeños pero con menor precisión.

---

## 7. Medidas del Error

Sea $A$ el **valor exacto** y $a$ la **aproximación**.

| Medida | Fórmula | Interpretación |
|--------|---------|----------------|
| Error absoluto | $\xi_a = \|A - a\|$ | Diferencia directa (en las unidades del problema) |
| Error relativo | $\delta_a = \xi_a / \|A\|$ | Proporción del error respecto al valor exacto |
| Error relativo porcentual | $\varepsilon_a = \delta_a \times 100\%$ | Lo mismo expresado en porcentaje |

**Ejemplo numérico:**
Si $A = 3.141592653...$ y $a = 3.1416$, entonces:
- $\xi_a = |3.141592... - 3.1416| = 0.0000073...$
- $\delta_a = 0.0000073.../3.141592... \approx 2.3 \times 10^{-6}$
- $\varepsilon_a \approx 0.00023\%$

---

## 8. Efectos del Redondeo en Operaciones Aritméticas

### Suma de números con magnitudes muy distintas
Al sumar $1 + 0.00001$:
```
(1)_{10}       = (0.1000 0000 0000 0000 0000 0000)_2 × 2^1
(0.00001)_{10} = (0.1010 0111 1100 0101 1010 1100)_2 × 2^{-16}
```
Resultado: $1.0000100136$ (hay pérdida de dígitos al alinear los exponentes)

### Resta de números cercanos (cancelación catastrófica)
La resta $1.00001 - 1$ genera pérdida de cifras significativas porque los dígitos más importantes se cancelan, amplificando el error relativo en el resultado.

**Ejemplo del límite de la derivada:**
$$d = \frac{\sin(1 + \theta) - \sin(1)}{\theta} \xrightarrow{\theta \to 0} \cos(1) \approx 0.54030$$

Al hacer $\theta$ muy pequeño en la computadora, el error aumenta en vez de disminuir — la resta de dos números casi iguales destruye las cifras significativas.

---

## 9. Resumen para la Exposición (10 minutos)

| Tiempo | Contenido |
|--------|-----------|
| 0–1 min | Presentación: ¿qué son los errores numéricos y por qué importan? |
| 1–3 min | Las 8 fuentes de error (énfasis en truncamiento y redondeo) |
| 3–5 min | Representación IEEE 754: bits, epsilon de máquina |
| 5–7 min | Overflow y underflow: explicación + resultados de tu máquina |
| 7–9 min | Error absoluto, relativo y porcentual |
| 9–10 min | Demo en vivo con R + conclusiones |

---

## 10. Preguntas que puede hacer el profesor

1. **¿Por qué el epsilon de máquina es aproximadamente $10^{-16}$?**
   Porque la mantisa tiene 52 bits, y $2^{-52} \approx 2.22 \times 10^{-16}$.

2. **¿Qué pasa si sumas el epsilon de máquina con 2 en vez de 1?**
   El intervalo entre números representables es mayor, así que necesitas un número mayor que $2\varepsilon$.

3. **¿Por qué el underflow tardó 1075 pasos y no 1022?**
   Por los números subnormales (denormalizados), que extienden el rango pero con menos precisión.

4. **¿Cómo afectan estos errores a un método como Newton-Raphson?**
   Pueden causar que el algoritmo nunca converja o que oscile si el error de redondeo domina sobre el criterio de parada.
