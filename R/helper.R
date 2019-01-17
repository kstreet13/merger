library(ggplot2)
library(dplyr)
library(tidyr)
plotARIs <- function(ARI) {
  p <- ARI %>% as.data.frame() %>%
          mutate(label = rownames(ARI)) %>%
          gather(key = label2, value = ari, -(ncol(ARI) + 1)) %>%
    ggplot(aes(x = label, y = label2, fill = ari)) +
    geom_tile() +
    guides(fill = F) +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(axis.line = element_blank())
  return(p)
}
