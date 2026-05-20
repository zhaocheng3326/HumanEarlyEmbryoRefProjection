# Example: project a query dataset with the package API.
# Run from the package source directory after installation:
#   R CMD INSTALL .
#   Rscript projecting_query_dataset.R

suppressPackageStartupMessages({
  library(HuEmProjection)
  library(dplyr)
  library(ggplot2)
})

predict_out <- project_demo_dataset(n = "n1000", run_milor = TRUE, cor_cutoff = 0.5)

print(
  predict_out$full.anno %>%
    group_by(pred_EML) %>%
    summarise(nCell = n_distinct(query_cell), .groups = "drop") %>%
    arrange(desc(nCell))
)

p <- ggplot() +
  geom_point(
    data = predict_out$umap %>% filter(pj != "query"),
    mapping = aes(x = UMAP_1, y = UMAP_2),
    color = "grey"
  ) +
  geom_point(
    data = predict_out$umap %>% filter(pj == "query"),
    mapping = aes(x = UMAP_1, y = UMAP_2),
    color = "red"
  )
print(p)

sessionInfo()
