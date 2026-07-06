library(tidyverse)
library(ggpubr)
library(readxl)

# 1. Load the data
df <- read_excel("data_fig_3g_colonic_organoids.xlsx")

# 2. Preparation: Rename "KO" to "Nr4a1 KO" so it matches your color scheme
df <- df %>%
  mutate(Genotype = ifelse(Genotype == "KO", "Nr4a1 KO", Genotype))

# 3. Ensure Genotype is a factor for consistent ordering (WT on the left)
df$Genotype <- factor(df$Genotype, levels = c("Nr4a1 KO", "WT"))

# 4. Create the plot
# Note: Since the column name has spaces and a %, we use backticks: `% of WT ave`
p_organoids <- ggplot(df, aes(x = Genotype, y = `% of WT ave`, color = Genotype, fill = Genotype)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, linewidth = 1, alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
  
  # Significance logic (Stars/ns)
  stat_compare_means(
    comparisons = list(c("WT", "Nr4a1 KO")),
    method = "t.test", 
    label = "p.signif",
    size = 6, 
    vjust = -0.5, 
    symnum.args = list(
      # Added 0.0001 as a threshold and the 4th asterisk
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  
  # Apply your specific colors
  scale_fill_manual(values = c("WT" = "#00BFC4", "Nr4a1 KO" = "#F8766D")) +
  scale_color_manual(values = c("WT" = "#00BFC4", "Nr4a1 KO" = "#F8766D")) +
  
  # Styling
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
  labs(
    #title = "Colonic Organoids",
    y = "% of WT average",
    x = "Genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), # Rotate x-axis labels
    axis.text.y = element_text(size = 14),
  )

# Display and save
print(p_organoids)
ggsave("plots/Colonic_Organoids_Plot.svg", p_organoids, width = 4, height = 4)
