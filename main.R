suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(CATALYST)
})

#tim::set_workflow_step_ids("https://tercen.com/tercen/w/bde20c9b529a415310cfcfa9b160c6fe/ds/3ac3de18-13e6-491b-b8b4-6c2a443b5a7c")

ctx = tercenCtx()

Normalization <- ctx$op.value("Normalization", as.character, "dvs")

data <- ctx$as.matrix() %>% t()
colnames(data) <- ctx$rselect()[[1]]

files <- ctx$cselect() %>% 
  select(contains("filename"))

fset <- data %>% as_tibble() %>% bind_cols(files) %>%
  group_by({if("filename" %in% names(.)) filename else NULL}) %>% 
  select(.,-filename)%>% 
  group_map(~tim::matrix_to_flowFrame(as.matrix(.x))) %>%
  flowCore::flowSet()

# construct SCE
sce <- prepData(fset)
# apply normalization; keep raw data
res <- normCytof(sce, beads =Normalization, k = 50, 
                 assays = c("counts", "exprs"), overwrite = TRUE, remove_beads = TRUE)
# check number & percentage of bead / removed events
# n <- ncol(sce); ns <- c(ncol(res$beads), ncol(res$removed))
# data.frame(
#   check.names = FALSE, 
#   "#" = c(ns[1], ns[2]), 
#   "%" = 100*c(ns[1]/n, ns[2]/n),
#   row.names = c("beads", "removed"))

# plot bead vs. dna scatters
sc_plot<-res$scatter
# plot smoothed bead intensities
line_plot<-res$lines 

sc_file <- suppressWarnings({tim::save_plot(sc_plot,
                                            type = "png",
                                            width = 750,
                                            height = 750,
                                            units = "px",
                                            dpi = 144,
                                            device = "png"
)})

line_file <- suppressWarnings({tim::save_plot(line_plot,
                                              type = "png",
                                              width = 750,
                                              height = 750,
                                              units = "px",
                                              dpi = 144,
                                              device = "png"
)})

#sc_df_plot <- tim::plot_file_to_df(sc_file) %>%
#  ctx$addNamespace() %>%
#  as_relation() 

#line_df_plot <- tim::plot_file_to_df(line_file) %>%
#  ctx$addNamespace() %>%
#  as_relation() 

df_plot<-bind_rows(tim::plot_file_to_df(sc_file),tim::plot_file_to_df(line_file))%>%
  ctx$addNamespace() %>%
  as_relation() 



# extract data excluding beads & doublets,
# and including normalized intensitied
sce <- res$data
#assayNames(sce)
#df <- assay(sce, "normexprs")
df <- assay(sce, "exprs")

rids <- ctx$rselect()[1]
colnames(rids) <- "channel"

df_out <- df %>%
  as_tibble(rownames = "channel") %>%
  tidyr::pivot_longer(cols = !contains("channel"), names_to = ".ci") %>%
  mutate(.ci = as.integer(gsub("V", "", .ci)) - 1L) %>%
  left_join(rids %>% mutate(.ri = seq(1, nrow(.)) - 1L), by = "channel") %>%
  ctx$addNamespace() %>%
  as_relation() 

join_res = df_out %>%
  left_join_relation(ctx$crelation, ".ci", ctx$crelation$rids) %>%
  left_join_relation(df_plot, list(), list()) %>%
  as_join_operator(ctx$cnames, ctx$cnames)

join_res %>%
  save_relation(ctx)

