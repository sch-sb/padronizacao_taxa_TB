install.packages("pacman")

install.packages("dsr")

pacman::p_load(tidyverse, dsr, janitor, readr)

sinan <- read_csv2("sinannet_cnv_tubercbr10092810_1_247_107.csv", 
                  locale = locale(encoding = "Latin1"))
pop <- 

dados <- sinan_tb %>%
  filter(ano == 2023) %>%
  mutate(faixa_etaria = cut(idade, breaks = c(0, 14, 44, 64, Inf),
                            labels = c("0-14", "15-44", "45-64", "65+")),
         sexo = case_when(sexo == "M" ~ "Masculino",
                          sexo == "F" ~ "Feminino")) %>%
  group_by(regiao_saude, faixa_etaria, sexo) %>%
  summarise(casos = n())


dados <- dados %>%
  left_join(populacao, by = c("regiao_saude", "faixa_etaria", "sexo"))

dados <- dados %>%
  left_join(pop_padrao, by = c("faixa_etaria", "sexo"))

resultados <- dsr::dsr(
  x = dados$casos,
  n = dados$pop,           # população observada
  stdpop = dados$pop_padrao,  # população padrão
  method = "gamma",
  group = dados$regiao_saude
)