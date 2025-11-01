# ============================================================
# Padronização de incidência (morbidade) de Tuberculose por idade
# Referência: https://epirhandbook.com/pt/new_pages/standardization.pt.html
# ============================================================

options(repos = c(CRAN = "https://cran.r-project.org"))

# -------------------------
# 0) Pacotes
# -------------------------
pkgs <- c("pacman")
if (!all(pkgs %in% rownames(installed.packages())))
  install.packages(setdiff(pkgs, rownames(installed.packages())))
suppressWarnings(suppressMessages(library(pacman)))

pacman::p_load(
  tidyverse, janitor, stringr, lubridate, readr, purrr, dsr, frailtypack, PHEindicatormethods)

# Instalação 
if (!requireNamespace("dsr", quietly = TRUE)) {
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/dsr/dsr_0.2.2.tar.gz"
  install.packages(packageurl, repos = NULL, type = "source")
}
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("PHEindicatormethods", quietly = TRUE)) {
  remotes::install_github("ukhsa-collaboration/PHEindicatormethods",
                          build_vignettes = FALSE, dependencies = TRUE)
}

options(scipen = 999)

# -------------------------
# 1) Parâmetros e caminhos
# -------------------------
PATH_BASE <- "C:/Users/sacha/Documents/Padronizacao_TB"
ARQ_TB    <- file.path(PATH_BASE, "dados_TB_tratado3.csv")
ARQ_POP   <- file.path(PATH_BASE, "pop_proj_is.csv")
ANO_INICIO <- 2013
ANO_FIM    <- 2024     
ANO_PADRAO <- 2022     

PATH_BASE <- "C:/Users/sacha/Documents/Padronizacao_TB"
ARQ_TB    <- file.path(PATH_BASE, "dados_TB_tratado3.csv")
ARQ_POP   <- file.path(PATH_BASE, "pop_proj_is.csv")

# -------------------------
# 1) Pacotes
# -------------------------
pkgs <- c("pacman")
if (!all(pkgs %in% rownames(installed.packages())))
  install.packages(setdiff(pkgs, rownames(installed.packages())))
suppressWarnings(suppressMessages(library(pacman)))

pacman::p_load(
  tidyverse,
  janitor,
  stringr,
  lubridate,
  readr,
  purrr,
  dsr,
  frailtypack,          
  PHEindicatormethods
)


if (!requireNamespace("dsr", quietly = TRUE)) {
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/dsr/dsr_0.2.2.tar.gz"
  install.packages(packageurl, repos = NULL, type = "source")
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("PHEindicatormethods", quietly = TRUE)) {
  remotes::install_github("ukhsa-collaboration/PHEindicatormethods",
                          build_vignettes = FALSE, dependencies = TRUE)
}

# -------------------------
# 2) Tabela de UFs 
# -------------------------
UF_NOMES <- tibble::tribble(
  ~cod_uf, ~sigla_uf, ~nome_uf,
  11L, "RO", "Rondônia",
  12L, "AC", "Acre",
  13L, "AM", "Amazonas",
  14L, "RR", "Roraima",
  15L, "PA", "Pará",
  16L, "AP", "Amapá",
  17L, "TO", "Tocantins",
  21L, "MA", "Maranhão",
  22L, "PI", "Piauí",
  23L, "CE", "Ceará",
  24L, "RN", "Rio Grande do Norte",
  25L, "PB", "Paraíba",
  26L, "PE", "Pernambuco",
  27L, "AL", "Alagoas",
  28L, "SE", "Sergipe",
  29L, "BA", "Bahia",
  31L, "MG", "Minas Gerais",
  32L, "ES", "Espírito Santo",
  33L, "RJ", "Rio de Janeiro",
  35L, "SP", "São Paulo",
  41L, "PR", "Paraná",
  42L, "SC", "Santa Catarina",
  43L, "RS", "Rio Grande do Sul",
  50L, "MS", "Mato Grosso do Sul",
  51L, "MT", "Mato Grosso",
  52L, "GO", "Goiás",
  53L, "DF", "Distrito Federal",
  99L, "BR", "Brasil"
)

# -------------------------
# 3) Funções auxiliares
# -------------------------


idade_to_gi <- function(idade_num) {
  gi <- floor(idade_num / 5) * 5
  gi <- ifelse(idade_num %in% 1:4, 1L, gi)
  gi[is.na(idade_num)] <- NA_integer_
  gi <- ifelse(gi > 80, 80L, gi)
  gi
}

norm_sexo_tb <- function(x) {
  x2 <- str_to_lower(x)
  case_when(
    str_detect(x2, "^f") ~ "Feminino",
    str_detect(x2, "^m") ~ "Masculino",
    TRUE ~ NA_character_
  )
}

norm_sexo_pop <- function(x) {
  x2 <- str_to_lower(str_trim(x))
  case_when(
    x2 %in% c("h", "homem", "masc", "masculino") ~ "Masculino",
    x2 %in% c("m", "mulher", "f", "fem", "feminino") ~ "Feminino",
    TRUE ~ NA_character_
  )
}

mun_to_uf <- function(cod_mun) {
  suppressWarnings(as.integer(floor(as.numeric(cod_mun) / 1000)))
}

# -------------------------
# 4) Leitura das bases
# -------------------------

# 4.1) TB (casos)
TB_raw <- read_csv2(ARQ_TB, show_col_types = FALSE) %>%
  clean_names()

TB <- TB_raw %>%
  mutate(
    ano        = coalesce(nu_ano, year(dt_notific), year(dt_noti_at), year(dt_diag)),
    sexo       = norm_sexo_tb(cs_sexo),
    cod_mun    = coalesce(id_mn_resi, id_municip, id_munic_a, id_munic_2),
    cod_uf     = coalesce(sg_uf, sg_uf_2, as.numeric(mun_to_uf(cod_mun))),
    idade_anos = coalesce(nu_idade_n_anos_completa, nu_idade_n_anos),
    gi         = idade_to_gi(idade_anos)
  ) %>%
  mutate(
    eh_caso_novo = case_when(
      "casos_novos" %in% names(.) ~ str_to_lower(casos_novos) == "sim",
      "tratamento"  %in% names(.) ~ str_detect(str_to_lower(tratamento), "caso\\s*novo"),
      TRUE ~ TRUE
    )
  ) %>%
  filter(!is.na(ano),
         ano >= ANO_INICIO,
         ano <= ANO_FIM,
         !is.na(gi), !is.na(sexo), !is.na(cod_uf),
         eh_caso_novo) %>%
  mutate(cod_uf = as.integer(cod_uf)) %>%
  dplyr::select(ano, cod_uf, sexo, gi)

# 4.2) População 
pop_raw <- read_csv2(ARQ_POP, show_col_types = FALSE,
                     locale = locale(encoding = "UTF-8")) %>%
  clean_names()

pop_uf <- pop_raw %>%
  mutate(
    sexo   = norm_sexo_pop(sexo),
    gi     = idade_to_gi(idade),
    cod_uf = as.integer(cod_uf),
    ano    = as.integer(ano)
  ) %>%
  filter(!is.na(ano),
         ano >= ANO_INICIO,
         ano <= ANO_FIM,
         !is.na(gi), !is.na(sexo), !is.na(cod_uf)) %>%
  group_by(ano, cod_uf, sexo, gi) %>%
  summarise(pop = sum(populacao, na.rm = TRUE), .groups = "drop")


POP_UF <- pop_uf
POP_BR <- POP_UF %>%
  group_by(ano, sexo, gi) %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
  mutate(cod_uf = 99L)

POP <- bind_rows(POP_UF, POP_BR) %>%
  arrange(ano, cod_uf, sexo, gi)

# -------------------------
# 5) Eventos agregados 
# -------------------------
CASOS_UF <- TB %>%
  group_by(ano, cod_uf, sexo, gi) %>%
  summarise(eventos = n(), .groups = "drop")

CASOS_BR <- TB %>%
  group_by(ano, sexo, gi) %>%
  summarise(eventos = n(), .groups = "drop") %>%
  mutate(cod_uf = 99L)

CASOS <- bind_rows(CASOS_UF, CASOS_BR) %>%
  arrange(ano, cod_uf, sexo, gi)

# -------------------------
# 6) Grade completa
# -------------------------
GRADE <- POP %>% distinct(ano, cod_uf, sexo, gi)

db_full <- GRADE %>%
  left_join(POP,  by = c("ano", "cod_uf", "sexo", "gi")) %>%
  left_join(CASOS, by = c("ano", "cod_uf", "sexo", "gi")) %>%
  mutate(eventos = replace_na(eventos, 0L))

# -------------------------
# 7) Dados - sexo
# -------------------------
db_ambos <- db_full %>%
  group_by(cod_uf, ano, gi) %>%
  summarise(
    populacao = sum(pop, na.rm = TRUE),
    eventos   = sum(eventos, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    sexo  = "ambos",
    grupo = paste0(stringr::str_pad(cod_uf, 2, pad = "0"), "-", ano)
  ) %>%
  dplyr::select(grupo, gi, sexo, populacao, eventos, cod_uf, ano)

# -------------------------
# 8) População padrão 
# -------------------------
if (!(ANO_PADRAO %in% unique(pop_raw$ano))) {
  ANO_PADRAO <- max(pop_raw$ano, na.rm = TRUE)
  message("ANO_PADRAO não encontrado na pop. Usando ano máximo disponível: ", ANO_PADRAO)
}

pop_padrao <- pop_raw %>%
  mutate(
    sexo = norm_sexo_pop(sexo),
    gi   = idade_to_gi(idade)
  ) %>%
  filter(ano == ANO_PADRAO) %>%
  group_by(gi) %>%
  summarise(pop = sum(populacao, na.rm = TRUE), .groups = "drop") %>%
  mutate(sexo = "ambos")

# -------------------------
# 9) Padronização com dsr
# -------------------------
TBI_padronizada_dsr <- dsr::dsr(
  data     = db_ambos,
  event    = eventos,
  fu       = populacao,
  subgroup = grupo,
  gi,
  sexo,
  refdata  = pop_padrao,
  method   = "gamma",
  sig      = 0.95,
  mp       = 100000,
  decimals = 2
)

# -------------------------
# 10) Limpeza dos nomes
# -------------------------
find_col <- function(df, ...) {
  pats <- c(...)
  nm <- names(df)
  idx <- Reduce(`&`, lapply(pats, function(p) grepl(p, nm, ignore.case = TRUE)))
  hit <- nm[idx]
  if (length(hit) == 0) {
    stop("Não achei coluna com padrões: ", paste(pats, collapse = " & "),
         "\nNomes disponíveis:\n", paste(nm, collapse = "\n"))
  }
  hit[1]
}

sub_col <- if ("Subgroup" %in% names(TBI_padronizada_dsr)) "Subgroup" else
  if ("subgroup" %in% names(TBI_padronizada_dsr)) "subgroup" else
    if ("Group" %in% names(TBI_padronizada_dsr)) "Group" else
      if ("group" %in% names(TBI_padronizada_dsr)) "group" else
        stop("Coluna de grupo não encontrada no resultado do dsr().")

TBI_padronizada_dsr <- TBI_padronizada_dsr %>%
  mutate(
    cod_uf = as.integer(substr(.data[[sub_col]], 1, 2)),
    ano    = readr::parse_integer(stringr::str_extract(.data[[sub_col]], "(\\d{4})$"))
  ) %>%
  # AQUI entra o teto de ano de novo
  filter(ano >= ANO_INICIO, ano <= ANO_FIM)

col_crude     <- find_col(TBI_padronizada_dsr, "crude", "rate")
col_crude_lcl <- find_col(TBI_padronizada_dsr, "(95|ci)", "(lc[il]|lower)", "crude")
col_crude_ucl <- find_col(TBI_padronizada_dsr, "(95|ci)", "(uc[il]|upper)", "crude")

col_std       <- find_col(TBI_padronizada_dsr, "(std|standard|age|adj|adjust)", "rate")
col_std_lcl   <- find_col(TBI_padronizada_dsr, "(95|ci)", "(lc[il]|lower)", "(std|standard|age|adj|adjust)")
col_std_ucl   <- find_col(TBI_padronizada_dsr, "(95|ci)", "(uc[il]|upper)", "(std|standard|age|adj|adjust)")

TBI_padronizada_dsr <- TBI_padronizada_dsr %>%
  mutate(
    TBI                 = .data[[col_crude]],
    TBI_limite_inferior = .data[[col_crude_lcl]],
    TBI_limite_superior = .data[[col_crude_ucl]],
    TBI_padronizada     = .data[[col_std]],
    TBI_padronizada_inf = .data[[col_std_lcl]],
    TBI_padronizada_sup = .data[[col_std_ucl]]
  ) %>%
  dplyr::select(
    cod_uf, ano,
    TBI, TBI_limite_inferior, TBI_limite_superior,
    TBI_padronizada, TBI_padronizada_inf, TBI_padronizada_sup
  ) %>%
  arrange(cod_uf, ano) %>%
  left_join(UF_NOMES, by = "cod_uf")

# -------------------------
# 11) Visualizações 
# -------------------------
uf_escolha <- 31L
dados_graf <- TBI_padronizada_dsr %>%
  filter(cod_uf == uf_escolha,
         !is.na(ano),
         ano >= ANO_INICIO,
         ano <= ANO_FIM)

titulo_uf <- if (nrow(dados_graf) > 0) {
  sprintf("Evolução da TBI e TBI padronizada para Tuberculose – %s (%s)",
          unique(dados_graf$nome_uf), unique(dados_graf$sigla_uf))
} else {
  sprintf("Evolução da TBI e TBI padronizada para Tuberculose – UF %02d", uf_escolha)
}

gg1 <- ggplot(dados_graf, aes(x = ano)) +
  geom_ribbon(aes(ymin = TBI_limite_inferior,
                  ymax = TBI_limite_superior,
                  fill = "TBI"), alpha = 0.25) +
  geom_ribbon(aes(ymin = TBI_padronizada_inf,
                  ymax = TBI_padronizada_sup,
                  fill = "TBI Padronizada"), alpha = 0.25) +
  geom_line(aes(y = TBI, color = "TBI")) +
  geom_line(aes(y = TBI_padronizada, color = "TBI Padronizada para Tuberculose")) +
  labs(
    title   = titulo_uf,
    caption = "Fonte: SINAN TB (casos novos) + Projeções por idade (IBGE).",
    y = "por 100.000 hab.", x = "Ano",
    fill = "Intervalos", color = "Taxas"
  ) +
  scale_x_continuous(
    limits = c(ANO_INICIO, ANO_FIM),
    breaks = seq(ANO_INICIO, ANO_FIM, by = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  ylim(c(0, NA)) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0))

print(gg1)

# -------------------------
# 11A) Gráficos por UF 
# -------------------------
OUT_DIR <- file.path(PATH_BASE, "graficos_uf")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

build_plot_uf <- function(uf_code) {
  dados_graf <- TBI_padronizada_dsr %>%
    filter(cod_uf == uf_code,
           !is.na(ano),
           ano >= ANO_INICIO,
           ano <= ANO_FIM)
  if (nrow(dados_graf) == 0) return(NULL)
  
  nome_uf  <- unique(dados_graf$nome_uf)
  sigla_uf <- unique(dados_graf$sigla_uf)
  titulo   <- sprintf("Evolução da TBI e TBI padronizada para Tuberculose – %s (%s)", nome_uf, sigla_uf)
  
  ggplot(dados_graf, aes(x = ano)) +
    geom_ribbon(aes(ymin = TBI_limite_inferior,
                    ymax = TBI_limite_superior,
                    fill = "TBI"), alpha = 0.25) +
    geom_ribbon(aes(ymin = TBI_padronizada_inf,
                    ymax = TBI_padronizada_sup,
                    fill = "TBI Padronizada"), alpha = 0.25) +
    geom_line(aes(y = TBI, color = "TBI")) +
    geom_line(aes(y = TBI_padronizada, color = "TBI Padronizada")) +
    labs(
      title   = titulo,
      caption = "Fonte: SINAN TB (casos novos) + Projeções por idade (IBGE).",
      y = "por 100.000 hab.", x = "Ano",
      fill = "Intervalos", color = "Taxas"
    ) +
    scale_x_continuous(
      limits = c(ANO_INICIO, ANO_FIM),
      breaks = seq(ANO_INICIO, ANO_FIM, by = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    ylim(c(0, NA)) +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0))
}

save_plot_uf <- function(uf_code) {
  g <- build_plot_uf(uf_code)
  if (is.null(g)) return(invisible(NULL))
  sigla <- UF_NOMES %>% filter(cod_uf == uf_code) %>% pull(sigla_uf)
  if (length(sigla) == 0 || is.na(sigla)) sigla <- sprintf("%02d", uf_code)
  arq <- file.path(OUT_DIR, sprintf("TBI_UF_%02d_%s.png", uf_code, sigla))
  ggsave(arq, g, width = 26, height = 14, units = "cm", dpi = 300, bg = "white")
  message(sprintf("Salvo: %s", arq))
  invisible(g)
}

ufs <- sort(unique(TBI_padronizada_dsr$cod_uf))
ufs <- setdiff(ufs, 99L)
walk(ufs, ~ try(save_plot_uf(.x), silent = TRUE))

# PDF único
pdf(file.path(OUT_DIR, "TBI_todas_ufs.pdf"), width = 11.69, height = 8.27)
walk(ufs, ~ {
  g <- suppressWarnings(try(build_plot_uf(.x), silent = TRUE))
  if (inherits(g, "ggplot")) print(g)
})
dev.off()

# -------------------------
# 12) Exportações
# -------------------------
write_csv(TBI_padronizada_dsr, file.path(PATH_BASE, "TBI_padronizada_dsr.csv"))
saveRDS(TBI_padronizada_dsr, file.path(PATH_BASE, "TBI_padronizada_dsr.rds"))
