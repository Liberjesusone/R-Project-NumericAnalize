# overflow
x <- 1.0
counter <- 0
while (is.finite(x)) {
    x <- x * 2
    counter <- counter + 1
}
cat("Desbordamiento ocurrio despues de", counter, "multiplicaciones\n")
cat("Maximo numero representable (.Machine):", .Machine$double.xmax, "\n")

# underflow
y <- 1.0
counter2 <- 0
while (y > 0) {
    y <- y / 2
    counter2 <- counter2 + 1
}
cat("Underflow ocurrio despues de", counter2, "divisiones\n")
cat("Minimo numero positivo (.Machine):", .Machine$double.xmin, "\n")
