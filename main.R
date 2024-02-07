suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(CATALYST)
})

ctx = tercenCtx()

beads <- ctx$op.value("beads", as.character, "dvs")
plot_width <- ctx$op.value("plot_width", as.double, 500)
plot_height <- ctx$op.value("plot_height", as.double, 500)

data <- ctx$as.matrix() %>% t()
row_data <- ctx$rselect()[[1]]
colnames(data) <- row_data

colnames(data) <- sub('^([0-9][0-9]+[0-9])([A-Z]+[a-z])', '\\2\\1', colnames(data))
colnames(data) <- sub('_.*$', '', colnames(data))

files <- ctx$cselect() %>% 
  select(contains("filename"))

# construct flowset
fset <- data %>% 
  as_tibble() %>% 
  bind_cols(files) %>%
  group_by({if("filename" %in% names(.)) filename else NULL}) %>% 
  select(.,-filename)%>% 
  group_map(~tim::matrix_to_flowFrame(as.matrix(.x))) %>%
  flowCore::flowSet()

# construct SCE
sce <- prepData(fset)
# apply normalization; replace raw data & remove beads
res <- normCytof(
  sce,
  beads = beads,
  k = 500, 
  assays = c("counts", "exprs"),
  overwrite = TRUE,
  remove_beads = FALSE
)

# plot bead vs. dna scatters
sc_file <- tim::save_plot(
  res$scatter,
  type = "png",
  width = plot_width,
  height = plot_height,
  units = "px",
  dpi = 144,
  device = "png"
)

# plot smoothed bead intensities
line_file <- tim::save_plot(
  res$lines,
  type = "png",
  width = plot_width,
  height = plot_height,
  units = "px",
  dpi = 144,
  device = "png"
)

#bind both plot
df_plot <- bind_rows(
    tercen::file_to_tercen(sc_file, filename = "Beads_vs_DNA_Scatters.png"),
    tercen::file_to_tercen(line_file, filename = "Smoothed_Bead_Intensities.png")
  ) %>%
  ctx$addNamespace() %>%
  as_relation(relation_name = "Diagnostic Plots")  %>%
  as_join_operator(list(), list())

# extract data excluding beads & doublets,
# and including normalized intensitied
df <- assay(res$data, "exprs")

rids <- seq_along(row_data) - 1L
names(rids) <- row_data

rownames(df) <- rids[rownames(df)]
colnames(df) <- as.integer(1:ncol(df) - 1)

df_out <- df %>%
  as_tibble(rownames = ".ri_norm") %>%
  tidyr::pivot_longer(cols = !contains(".ri_norm"), names_to = ".ci_norm") %>%
  mutate(.ri_norm = as.integer(.ri_norm), .ci_norm = as.integer(.ci_norm)) %>%
  ctx$addNamespace() %>%
  as_relation(relation_name = "Normalised Data")

df_beads <- tibble(is_bead = as.double(res$data$is_bead)) %>%
  mutate(.ci_beads = 1:nrow(.) - 1L) %>%
  ctx$addNamespace() %>%
  as_relation(relation_name = "Beads")

join_res = df_out %>%
  left_join_relation(df_beads, ".ci_norm", ".ci_beads") %>%
  left_join_relation(ctx$rrelation, ".ri_norm", ctx$rrelation$rids) %>%
  left_join_relation(ctx$crelation, ".ci_norm", ctx$crelation$rids) %>%
  as_join_operator(c(ctx$rnames, ctx$cnames), c(ctx$rnames, ctx$cnames))

result <- save_relation(list(join_res, df_plot), ctx)

