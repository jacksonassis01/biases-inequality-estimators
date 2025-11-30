# Arquivo que inicializa testthat
library(testthat)

# Carrega todos os arquivos da pasta R/
for (f in list.files("../R", full.names = TRUE)) {
  source(f)
}

test_check("main")
