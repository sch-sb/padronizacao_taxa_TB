install.packages("pacman")


pacman::p_load(tidyverse, dsr, janitor, readr, readxl)

pop_raw <- read_xlsx("projecoes_2024_tab1_idade_simples.xlsx",
                    skip = 4) %>% clean_names() 

pop_2024 <- pop_raw %>%
  select(x1, x2, x3, x4, x5, x30) %>%
  filter(!x5 %in% c("Brasil", "Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste")) %>%
  filter(x2 == "Ambos") %>%
  mutate(x1 = as.numeric(x1)) %>%
  filter(!is.na(x1))


pop_agrupada <- pop_2024 %>%
  rename(
    idade = x1,
    sexo = x2,
    cod_uf = x3,
    sigla_uf = x4,
    uf = x5,
    pop = x30
  ) %>%
  mutate(faixa_etaria = cut(idade,
                            breaks = c(-1, 14, 44, 64, Inf),
                            labels = c("0-14", "15-44", "45-64", "65+"))) %>%
  group_by(uf, faixa_etaria) %>%
  summarise(pop = sum(pop), .groups = 'drop')



sinan_raw <- read_csv2("sinannet_cnv_tubercbr10092810_1_247_107.csv",
                       locale = locale(encoding = "Latin1")) %>% clean_names()

ufs_map <- data.frame(
  uf_cod = c("11", "12", "13", "14", "15", "16", "17", "21", "22", "23", "24", "25", "26", "27", "28", "29", "31", "32", "33", "35", "41", "42", "43", "50", "51", "52", "53"),
  uf_nome = c("Rondônia", "Acre", "Amazonas", "Roraima", "Pará", "Amapá", "Tocantins", "Maranhão", "Piauí", "Ceará", "Rio Grande do Norte", "Paraíba", "Pernambuco", "Alagoas", "Sergipe", "Bahia", "Minas Gerais", "Espírito Santo", "Rio de Janeiro", "São Paulo", "Paraná", "Santa Catarina", "Rio Grande do Sul", "Mato Grosso do Sul", "Mato Grosso", "Goiás", "Distrito Federal")
)

sinan_com_uf <- sinan_raw %>%
  mutate(uf_cod = str_sub(regiao_de_saude_cir_de_notif, 1, 2)) %>%
  left_join(ufs_map, by = "uf_cod")

dados_casos <- sinan_com_uf %>%
  # Seleciona as colunas relevantes
  select(uf = uf_nome, starts_with("x")) %>%
  pivot_longer(cols = -uf,
               names_to = "grupo_etario_original",
               values_to = "casos") %>%
  mutate(casos = as.numeric(replace(casos, casos == "-", NA)),
         casos = coalesce(casos, 0)) %>%
  mutate(faixa_etaria = case_when(
    grupo_etario_original %in% c("x_1_ano", "x01_04", "x05_09", "x10_14") ~ "0-14",
    grupo_etario_original %in% c("x15_19", "x20_39") ~ "15-44",
    grupo_etario_original %in% c("x40_59", "x60_64") ~ "45-64",
    grupo_etario_original %in% c("x65_69", "x70_79", "x80_e") ~ "65+"
  )) %>%
  filter(!is.na(faixa_etaria) & !is.na(uf)) %>%
  group_by(uf, faixa_etaria) %>%
  summarise(casos = sum(casos), .groups = 'drop')


dados_completos <- left_join(dados_casos, pop_agrupada, by = c("uf", "faixa_etaria"))


pop_padrao <- pop_agrupada %>%
  group_by(faixa_etaria) %>%
  summarise(pop_padrao = sum(pop))


dados_finais <- left_join(dados_completos, pop_padrao, by = "faixa_etaria") %>%
  filter(!is.na(pop) & !is.na(casos)) 


resultados_normalizados <- dsr::dsr(
  data = dados_finais,
  event = casos,
  pop = pop,
  stdpop = pop_padrao,
  group = uf,
  method = "gamma",
  ci.level = 0.95
)


print(resultados_normalizados)
