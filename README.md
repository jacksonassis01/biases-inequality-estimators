# Vieses em estimadores de Desigualdade

Este projeto implementa simulações Monte Carlo para análise de índices de desigualdade (Theil, Atkinson e o índice de dispersão (VMR) em populações geradas a partir de misturas de distribuições Gama com parâmetro de taxa constante. Também inclui aplicação em dados reais de PIB dos países norte americanos e oceânicos.

------------------------------------------------------------------------

## Estrutura do Projeto

```         
.
├── R/                          # Funções do projeto
│   ├── indices/                # Cálculo de índices e expectativas teóricas
│   ├── mixture/                # Geração e ajuste de misturas
│   └── montecarlo/             # Setup e execução do experimento Monte Carlo
├── application-real-data/      # Aplicação dos métodos em dados reais
│   ├── data/
│   └── figures/
├── dissertation/               # Scripts para gerar figuras e tabelas da dissertação
├── results/                    # Resultados salvos de simulações
├── tests/                      # Testes unitários
└── main.R                      # Script principal para rodar o projeto
```

------------------------------------------------------------------------

## Dependências

O projeto utiliza os seguintes pacotes R (instale caso não tenha):

``` r
install.packages(c(
  "dplyr",
  "ggplot2",
  "scales",
  "reshape2",
  "patchwork",
  "stringr",
  "hypergeo",
  "testthat"
))
```

> Se você usar `tidyverse`, ele inclui `dplyr`, `ggplot2` e `readr`.

------------------------------------------------------------------------

## Como Rodar

1.  Abra o RStudio ou R no diretório raiz do projeto.
2.  Carregue todas as funções e rode a simulação principal:

``` r
source("main.R")
```

O script `main.R` faz o seguinte:

-   Carrega todas as funções do projeto (`R/indices`, `R/mixture`, `R/montecarlo`).
-   Roda a simulação Monte Carlo com os parâmetros definidos.
-   Calcula as métricas da simulação (Média, EQM, Erro Padrão (MCSE)).
-   Salva resultados em `results/`.
-   Pode gerar figuras e tabelas na pasta `dissertation/`.
-   Pode aplicar os métodos em dados reais (`application-real-data/`).

------------------------------------------------------------------------

## Configuração de Parâmetros

No `main.R` você pode alterar:

-   `pii` → probabilidades da mistura.
-   `alpha1_values`, `alpha2_values` → parâmetros dos índices.
-   `beta` → parâmetro de escala da mistura.
-   `n_values` → tamanhos de amostra.
-   `Nrep` → número de repetições Monte Carlo.

------------------------------------------------------------------------

## Resultados

Os resultados são salvos em `.rds` na pasta `results/`. Por exemplo:

-   `R_1000_n_10_100_alpha1_0.50_alpha2_0.50_5_beta_1.00.rds` contém os resultados do experimento com 1000 repetições.
-   `results_with_metrics.rds` inclui métricas MCSE.

Para carregar:

``` r
res <- readRDS("results/results_with_metrics.rds")
```

------------------------------------------------------------------------

## Dissertação

Os scripts em `dissertation/` geram:

-   Figuras de viés e MCSE.
-   Tabelas de MCSE.

Use:

``` r
source("dissertation/main.R")
```

------------------------------------------------------------------------

## Aplicação em Dados Reais

Os dados estão em `application-real-data/data/`. O script `application-real-data/main.R` processa os dados e gera figuras em `application-real-data/figures/`.

------------------------------------------------------------------------

## Testes

Para rodar testes unitários:

``` r
library(testthat)
test_dir("tests/testthat")
```

------------------------------------------------------------------------

## Observações

-   Recomenda-se usar o RStudio para facilitar a execução.
-   Os arquivos `.rds` permitem salvar e carregar resultados sem precisar refazer simulações longas.
-   Certifique-se de ter todas as dependências instaladas antes de rodar o `main.R`.

------------------------------------------------------------------------
