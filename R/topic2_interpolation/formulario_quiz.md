# Formulario — Interpolación Polinomial

---

## Lagrange

$$P(x) = \sum_{k=0}^{n} L_k(x)\, f(x_k)$$

$$L_k(x) = \prod_{\substack{j=0\\j\neq k}}^{n} \frac{x - x_j}{x_k - x_j}$$

---

## Newton hacia adelante · datos equiespaciados, evaluar cerca de $x_0$

$$\Delta^0 f_i = f_i \qquad \Delta^k f_i = \Delta^{k-1}f_{i+1} - \Delta^{k-1}f_i$$

$$s = \frac{x - x_0}{h} \qquad P(x) = \sum_{k=0}^{n} \binom{s}{k} \Delta^k f_0$$

$$\binom{s}{k} = \frac{s(s-1)(s-2)\cdots(s-k+1)}{k!}$$

Forma explícita:

$$P(x) = f_0 + s\,\Delta f_0 + \frac{s(s-1)}{2!}\Delta^2 f_0 + \frac{s(s-1)(s-2)}{3!}\Delta^3 f_0 + \cdots$$

---

## Newton hacia atrás · datos equiespaciados, evaluar cerca de $x_n$

$$\nabla^0 f_i = f_i \qquad \nabla^k f_i = \nabla^{k-1}f_i - \nabla^{k-1}f_{i-1}$$

$$s = \frac{x - x_n}{h} \qquad P(x) = f_n + s\,\nabla f_n + \frac{s(s+1)}{2!}\nabla^2 f_n + \frac{s(s+1)(s+2)}{3!}\nabla^3 f_n + \cdots$$

---

## Newton diferencias divididas · cualquier espaciado

$$f[x_i] = f(x_i) \qquad f[x_i,\ldots,x_{i+k}] = \frac{f[x_{i+1},\ldots,x_{i+k}] - f[x_i,\ldots,x_{i+k-1}]}{x_{i+k} - x_i}$$

$$P(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \cdots$$

Coeficientes = **primera fila** de la tabla triangular.

---

## Error de truncamiento

$$E(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^{n}(x - x_i) \qquad \xi \in \bigl[\min_i x_i,\; \max_i x_i\bigr]$$

---

## Hermite · usa $f(x_k)$ **y** $f'(x_k)$

Nodos dobles: $z_{2k} = z_{2k+1} = x_k$. Diferencia dividida coincidente:

$$f[z_{2k}, z_{2k+1}] = f'(x_k)$$

Luego aplicar **Newton diferencias divididas** sobre los $z_i$.

$$E(x) = \frac{f^{(2n+2)}(\xi)}{(2n+2)!} \prod_{i=0}^{n}(x - x_i)^2$$

---

## Chebyshev · nodos óptimos en $[a,\, b]$

$$x_k = \frac{a+b}{2} + \frac{b-a}{2}\cos\!\left(\frac{2k+1}{2(n+1)}\,\pi\right), \quad k = 0, 1, \ldots, n$$

Cota mínima del producto:

$$\max_{x \in [a,b]}\left|\prod_{i=0}^{n}(x-x_i)\right| = \frac{(b-a)^{n+1}}{2^{2n+1}}$$

---

## Splines cúbicos naturales

$h_i = x_{i+1} - x_i$. En cada intervalo $[x_i, x_{i+1}]$:

$$S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3$$

$$a_i = f_i \qquad c_i = \frac{M_i}{2} \qquad d_i = \frac{M_{i+1}-M_i}{6h_i} \qquad b_i = \frac{f_{i+1}-f_i}{h_i} - \frac{h_i(2M_i+M_{i+1})}{6}$$

Sistema tridiagonal para $M_i = S''(x_i)$:

$$h_{i-1}M_{i-1} + 2(h_{i-1}+h_i)M_i + h_iM_{i+1} = 6\!\left(\frac{f_{i+1}-f_i}{h_i} - \frac{f_i-f_{i-1}}{h_{i-1}}\right)$$

Condiciones de frontera naturales: $M_0 = 0,\quad M_n = 0$
