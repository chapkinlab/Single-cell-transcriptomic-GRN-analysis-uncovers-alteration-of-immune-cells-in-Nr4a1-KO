library(tidyverse)
library(ggpubr)
library(patchwork)

# 1. Prepare Data
df <- read.csv("nr4a1_qpcr.csv") %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "KO")),
    Tissue   = factor(Tissue, levels = c("Colon", "Liver")),
    Gene     = factor(Gene)
  )

# 2. Refined Plotting Function
plot_qpcr <- function(data, gene_name) {
  
  plot_data <- data %>% filter(Gene == gene_name)
  
  ggplot(plot_data, aes(x = Genotype, y = Expression, color = Genotype, fill = Genotype)) +
    stat_summary(fun = mean, geom = "bar", width = 0.6, linewidth = 1, alpha = 0.2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
    geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
    facet_wrap(~ Tissue, scales = "free_y") +
    
    # Re-adding the significance logic
    stat_compare_means(
      comparisons = list(c("WT", "KO")),
      method = "t.test", 
      label = "p.signif",
      vjust = -0.5, # Adjusted to ensure stars are visible above the bracket
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
        symbols = c("***", "**", "*", "ns")
      )
    ) +
    
    scale_fill_manual(values = c("WT" = "#4C72B0", "KO" = "#8172B3")) +
    scale_color_manual(values = c("WT" = "#4C72B0", "KO" = "#8172B3")) +
    
    # Expand Y axis so stars don't get cut off at the top
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    
    labs(
      title = paste(gene_name, "expression"),
      y = "mRNA expression (normalized to TBP)",
      x = "Genotype"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold.italic", size = 16),
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 18)
    )
}

# 3. Generate and Label (A-D)
p_nr4a1 <- plot_qpcr(df, "Nr4a1")
p_nr4a2 <- plot_qpcr(df, "Nr4a2")

# Use '+' instead of '&' for the annotation to avoid overriding the stat_compare_means
final_figure <- (p_nr4a1 / p_nr4a2) + 
  plot_annotation(tag_levels = 'A')

# Apply the tag theme separately to keep it clean
final_figure <- final_figure & 
  theme(plot.tag = element_text(size = 22, face = "bold"))

final_figure
# 4. Save
ggsave("plots/Figure_S1_Final_Annotated.svg", final_figure, width = 7, height = 9)
