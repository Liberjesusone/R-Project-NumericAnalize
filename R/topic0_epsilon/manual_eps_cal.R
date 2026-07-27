eps <- 1.0
while ((eps / 2 ) + 1 > 1) {
  eps <- eps / 2
}
cat("Epsilon de maquina calculado:", eps, "\n")
cat("Epsilon de maquina en R (.Machine):", .Machine$double.eps, "\n")
