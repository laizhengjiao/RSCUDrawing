library(ggplot2)
library(ggpubr)
library(dplyr)

file_path <- "your_path/your_species_RSCU_stack.csv" 
one <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

one <- one %>%
  mutate(
    AA = case_when(
      AA %in% c("Leu1", "Leu2") ~ "Leu",
      AA %in% c("Ser1", "Ser2") ~ "Ser",
      TRUE ~ AA
    )
  )

aa_freq <- one %>%
  filter(!is.na(aaRatio)) %>%
  group_by(AA) %>%
  summarise(aaRatio = sum(aaRatio), .groups = "drop")

one <- one %>%
  group_by(AA) %>%
  mutate(Fill = row_number(),
         CodonIndex = row_number()) %>%
  ungroup()

aa_order <- c("Gln","His","Asn","Pro","Thr",
              "Leu","Glu","Met","Arg","Tyr","Asp","Lys",
              "Ala","Ile","Ser",
              "Cys","Trp","Val","Gly","Phe")

fill_colors <- c("#4E79A7","#E05C4B","#76B041","#7B5EA7","#F0B429","#3AAFA9","#C96A65","#8BAE5A")
one$AA_factor <- factor(one$AA, levels = aa_order)
xfreq <- factor(aa_freq$AA, levels = aa_order)

base_theme <- theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
  axis.line.y = element_line(color = "black", linewidth = 0.5),
  axis.line.x = element_line(color = "black", linewidth = 0.5),
  axis.ticks = element_line(color = "black", linewidth = 0.3),
  axis.ticks.length = unit(3, "pt"),
  axis.text = element_text(size = 8, color = "black"),
  axis.title.y = element_text(size = 9, margin = margin(r = 5)),
  legend.position = "none"
)

p_freq <- ggplot(aa_freq, aes(x = xfreq, y = aaRatio)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#5B8DB8") +
  geom_text(aes(label = sprintf("%.2f", aaRatio)), vjust = -0.4, size = 2.1, color = "grey35") +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0,16,4), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "a", x = NULL, y = "Amino acid frequency (%)") +
  base_theme +
  theme(axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 19, hjust = 0, margin = margin(b = 4, l = -14)),
        plot.margin = margin(10, 10, 6, 10))

p_rscu <- ggplot(one, aes(x = AA_factor, y = RSCU)) +
  geom_col(aes(fill = as.factor(Fill)), width = 0.7, position = "stack") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(0, 9.5), breaks = seq(0,8,2), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "b", x = NULL, y = "RSCU") +
  base_theme +
  theme(axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 19, hjust = 0, margin = margin(b = 4, l = -14)),
        plot.margin = margin(0, 10, 0, 10))

cell_h <- 1
p_codon <- ggplot() +
  geom_rect(data = one,
            aes(xmin = as.numeric(AA_factor)-0.46, xmax = as.numeric(AA_factor)+0.46,
                ymin = -(CodonIndex-1)*cell_h, ymax = -CodonIndex*cell_h),
            fill = fill_colors[one$Fill], color = "black", linewidth = 0.3) +
  geom_text(data = one,
            aes(x = as.numeric(AA_factor), y = -(CodonIndex-0.5)*cell_h, label = Codon),
            size = 1.9, colour = "white", fontface = "bold") +
  scale_x_continuous(limits = c(0.5, length(aa_order)+0.5), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_void() +
  theme(plot.margin = margin(12,10,8,10))

p_b <- ggarrange(p_rscu, p_codon, heights = c(1.6, 0.45), ncol = 1, align = "v")
pall <- ggarrange(p_freq, p_b, heights = c(1, 1.4), ncol = 1)

ggsave("RSCU_AA_fixed.png", plot = pall, width = 10, height = 10, dpi = 900, bg = "white")
