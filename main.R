suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(CATALYST)
})

ctx = tercenCtx()

Normalization <- ctx$op.value("Normalization", as.character, "dvs")

data <- ctx$as.matrix() %>% t()
colnames(data) <- ctx$rselect()[[1]]

files <- ctx$cselect() %>% 
  select(contains("filename"))

# construct flowset
fset <- data %>% as_tibble() %>% bind_cols(files) %>%
  group_by({if("filename" %in% names(.)) filename else NULL}) %>% 
  select(.,-filename)%>% 
  group_map(~tim::matrix_to_flowFrame(as.matrix(.x))) %>%
  flowCore::flowSet()

# construct SCE
sce <- prepData(fset)
# apply normalization; replace raw data & remove beads
res <- normCytof(sce, beads =Normalization, k = 50, 
                 assays = c("counts", "exprs"), overwrite = TRUE, remove_beads = TRUE)

# plot bead vs. dna scatters
sc_plot<-res$scatter
sc_file <- suppressWarnings({tim::save_plot(sc_plot,
                                            type = "png",
                                            width = 750,
                                            height = 750,
                                            units = "px",
                                            dpi = 144,
                                            device = "png"
)})

# plot smoothed bead intensities
line_plot<-res$lines 
line_file <- suppressWarnings({tim::save_plot(line_plot,
                                              type = "png",
                                              width = 750,
                                              height = 750,
                                              units = "px",
                                              dpi = 144,
                                              device = "png"
)})

#bind both plot
df_plot<-bind_rows(tim::plot_file_to_df(sc_file),tim::plot_file_to_df(line_file))%>%
  ctx$addNamespace() %>%
  as_relation() 



# extract data excluding beads & doublets,
# and including normalized intensitied
sce <- res$data
df <- assay(sce, "exprs")

rids <- ctx$rselect()[1]
colnames(rids) <- "variable"

df_out <- df %>%
  as_tibble(rownames = "variable") %>%
  tidyr::pivot_longer(cols = !contains("variable"), names_to = ".ci") %>%
  mutate(.ci = as.integer(gsub("V", "", .ci)) - 1L) %>%
  left_join(rids %>% mutate(.ri = seq(1, nrow(.)) - 1L), by = "variable") %>%
  ctx$addNamespace() %>%
  as_relation() 

join_res = df_out %>%
  left_join_relation(ctx$crelation, ".ci", ctx$crelation$rids) %>%
  left_join_relation(df_plot, list(), list()) %>%
  as_join_operator(ctx$cnames, ctx$cnames)

join_res %>%
  save_relation(ctx)

