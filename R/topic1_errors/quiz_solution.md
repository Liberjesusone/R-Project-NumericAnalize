# Solución del Quiz de Errores en Métodos Numéricos
## Problema 1
Evalúe $f(x) = e^x - x$ en $x = 0.1$ usando serie de Taylor de tercer orden centrada en $x_0 = 0$ y $x_0 = 1$. Explique la diferencia.

### Derivadas de $f(x)$
$$f(x) = e^x - x \qquad f(x_0) \quad ; \quad f'(x) = e^x - 1 \qquad f'(x_0)$$
$$f''(x) = e^x \qquad f''(x_0) \quad ; \quad f'''(x) = e^x \qquad f'''(x_0)$$
La fórmula general de Taylor de tercer orden es: $$P_3(x) = f(x_0) + f'(x_0)(x - x_0) + \frac{f''(x_0)}{2!}(x-x_0)^2 + \frac{f'''(x_0)}{3!}(x-x_0)^3$$
### Centrado en $x_0 = 0$
- Evaluando en $x_0 = 0$: $$f(0) = 1, \quad f'(0) = 0, \quad f''(0) = 1, \quad f'''(0) = 1$$$$P_3(x) = 1 + 0 \cdot x + \frac{x^2}{2} + \frac{x^3}{6}$$
- Evaluando en $x = 0.1$: $$P_3(0.1) = 1 + \frac{(0.1)^2}{2} + \frac{(0.1)^3}{6} = 1 + 0.005 + 0.000167 \approx 1.005167$$ 
### Centrado en $x_0 = 1$
- Evaluando en $x_0 = 1$ (usando $e \approx 2.71828$): $$f(1) = e - 1 \approx 1.71828, \quad f'(1) = e - 1 \approx 1.71828, \quad f''(1) = e \approx 2.71828, \quad f'''(1) = e \approx 2.71828$$ 
- Sea $h = x - x_0 = 0.1 - 1 = -0.9$: $$P_3(0.1) = (e-1) + (e-1)(-0.9) + \frac{e\,(0.81)}{2} + \frac{e\,(-0.729)}{6}$$ $$= 1.71828 - 1.54645 + 1.10020 - 0.33048 \approx 0.9416$$
### Comparación con el valor exacto
$$f(0.1) = e^{0.1} - 0.1 = 1.10517 - 0.1 = 1.00517$$

| Centro    | $P_3(0.1)$ | Error absoluto               |
| --------- | ---------- | ---------------------------- |
| $x_0 = 0$ | $1.005167$ | $\approx 4 \times 10^{-6}$   |
| $x_0 = 1$ | $0.9416$   | $\approx 6.4 \times 10^{-2}$ |

### Explicación de la diferencia
La serie converge más rápido cuanto más cerca está $x$ del centro $x_0$. Como $x = 0.1$ está mucho más cerca de $0$ que de $1$, el polinomio centrado en $x_0 = 0$ es más preciso con el mismo número de términos.

### Diferencias con mi respuestas. 
- En mi quiz me equivoqué al calcular la primera derivada para $e^{ x } -x$ en por ende eso afecto un poco la obtención más precisa de los resultados 


---
## Problema 2
Procesador de 8 bits: 1 bit de signo, 3 bits de exponente, 4 bits de mantisa. Calcule el máximo y mínimo valor representable.

### Distribución de bits
```
 7  |  6  5  4  |  3  2  1  0
 S  |  E  E  E  |  M  M  M  M
```
- **S**: bit de signo (0 = positivo, 1 = negativo)
- **EEE**: exponente en notación sesgada (biased)
- **MMMM**: mantisa (parte fraccionaria; el 1 inicial es implícito en números normales)

### Sesgo del exponente (Bias)
Con $k = 3$ bits de exponente: $$\text{Bias} = 2^{k-1} - 1 = 2^2 - 1 = 3$$Los valores reservados son $000_2$ (cero y subnormales) y $111_2$ (Inf y NaN). La tabla de exponentes usables es:

| Almacenado | Decimal | Exponente real |
|------------|---------|----------------|
| $000$ | 0 | — (subnormal/cero) |
| $001$ | 1 | $1 - 3 = -2$ |
| $010$ | 2 | $2 - 3 = -1$ |
| $011$ | 3 | $3 - 3 = \phantom{-}0$ |
| $100$ | 4 | $4 - 3 = \phantom{-}1$ |
| $101$ | 5 | $5 - 3 = \phantom{-}2$ |
| $110$ | 6 | $6 - 3 = \phantom{-}3$ |
| $111$ | 7 | — (Inf/NaN) |

Rango de exponentes normales: $-2$ a $+3$.

### Valor máximo representable
Signo positivo, exponente máximo ($110_2 \to$ real $= 3$), mantisa máxima ($1111_2$): $$1.1111_2 = 1 + \frac{1}{2} + \frac{1}{4} + \frac{1}{8} + \frac{1}{16} = \frac{31}{16} = 1.9375$$
$$V_{\max} = 1.9375 \times 2^3 = 15.5$$
El mínimo negativo es $-15.5$ (mismo valor, bit de signo $= 1$).

### Mínimo positivo normal
Exponente mínimo normal ($001_2 \to$ real $= -2$), mantisa mínima ($0000_2$): $$V_{\min,\text{normal}} = 1.0000_2 \times 2^{-2} = 1 \times \frac{1}{4} = 0.25$$
### Mínimo positivo subnormal
Cuando el exponente almacenado es $000_2$, el exponente real se fija en $-2$ pero no hay 1 implícito — la mantisa empieza con $0$: $$0.0001_2 \times 2^{-2} = 2^{-4} \times 2^{-2} = 2^{-6}$$$$V_{\min,\text{sub}} = \frac{1}{64} = 0.015625$$
Los números subnormales permiten representar valores más pequeños que el mínimo normal a costa de perder bits de precisión progresivamente. Partimos del los bits 0 001 0000, donde el exponente es $001 = 1_{10}$ -> $1-3 = -2$, lo cuál es el número $1.000 * 2^{-2} = 0.25$, desde allí podemos dividir entre dos y básicamente se irán desplazando los bits hacia la derecha 
- $0~ 000 ~1000 = 0.125_{10}$
- $0~ 000 ~0100 = 0.0625._{10}$
- $0 ~000 ~0010 = 0.03125_{10}$
- $0 ~000 ~0001 = 0.015625_{10}$
- $0~ 000~ 0000 = 0_{10}$

### Resumen

| Cantidad                  | Valor      |
| ------------------------- | ---------- |
| Máximo positivo           | $15.5$     |
| Mínimo negativo           | $-15.5$    |
| Mínimo normal positivo    | $0.25$     |
| Mínimo subnormal positivo | $0.015625$ |

### Diferencias con mis respuestas 
- Al obtener el rango del exponente no tome en cuenta que eliminamos el primer 000 y el último 111 y que eso al dejaros solo 6 unidades, la resta para BIASED debía haber sido -3 y no -4 cómo yo elegí
- Con respecto al desplazamiento de bits para obtener los números subnormales, no usé la representación con el bit más pequeño diferente de 000, del cuál empezábamos a desplazar para obtener así el número más pequeño representable. 
- En la representación de la mantisa, el rango es de 0-15 cómo mencioné, pero ese no es un dígito que directamente se coloca tras la coma, cada bit de la mantisa, representa un multiplo de una potencia negativa de 2 