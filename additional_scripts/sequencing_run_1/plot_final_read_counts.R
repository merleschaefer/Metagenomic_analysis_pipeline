# ── Libraries ────────────────────────────────────────────
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
# ── Paths ────────────────────────────────────────────────
READS <- "/path/to/final_read_counts.tsv"
USI_DESCRIPTION <- "/path/to/USI_description.csv" #for translation of lab internal codes to sample descriptions for interpretability
VIRUS_ACCESSION_NAME <- "/path/to/contig_name.csv" #for translation of accession codes to virus names for plotting
OUTPUT_DIR <- "/path/to/output_directory"

# ── Load reads ────────────────────────────────────
reads <- read.table(file=READS, sep = '\t',header = TRUE)

# ── Shorten read names────────────────────────────────────
reads$sample <- sub("_.*", "", reads$sample)

# ── Load USI description and merge ────────────────────────────────────
usi_description <- read.table(file=USI_DESCRIPTION, sep = ',',header = TRUE)
reads <- reads %>%
  inner_join(usi_description, by=c("sample"="USI"))

# ── Load virus accession name and merge────────────────────────────────────
virus_accession_name <- read.table(file=VIRUS_ACCESSION_NAME, sep = ',',header = TRUE)
reads <- reads%>%
  inner_join(virus_accession_name, by=c("accession"="identifier"))


#── only keep viruses with reads ────────────────────────────────────
reads <- reads %>%
  group_by(species)%>%
  filter(sum(reads_without_dups)>0)%>%
  ungroup()


#── make sums of anelloviruses ────────────────────────────────────
prefixes <- c("Alphatorquevirus", "Betatorquevirus", "Gammatorquevirus")

pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")

reads_sum <- reads %>%
  mutate(
    matched_prefix = sub(paste0(pattern, ".*"), "\\1", species),
    is_match       = grepl(pattern, species),
    group_key      = ifelse(is_match, matched_prefix, as.character(row_number())),
    species_new    = ifelse(is_match, matched_prefix, species)  
  ) %>%
  group_by(sample, group_key) %>%
  summarise(
    species            = first(species_new),
    reads_with_dups    = sum(reads_with_dups),
    reads_without_dups = sum(reads_without_dups),
    across(-c(species_new, reads_with_dups, reads_without_dups, matched_prefix, is_match), first),
    .groups = "drop"
  ) %>%
  select(-group_key)


#── order viruses (by prevalence) and samples ────────────────────────────────────
order_viruses <- function(df) {
  virus_order <- df %>%
    group_by(species) %>%
    summarise(n_detected = sum(reads_without_dups>0)) %>%
    arrange(desc(n_detected)) %>%
    pull(species)
  
  virus_levels <- rev(virus_order)
  
  df <- df %>%
    mutate(species = factor(species, levels = virus_levels))
}

reads_sum  <- order_viruses(reads_sum)


SAMPLE_ORDER <- c(
  "Pos-V0-CF",
  "Pos-V2-CF", 
  "Pos-V20-CF",
  "Pos-V200-CF",
  "Pos-V8-CT",
  "Neg-V20-CF", 
  "Neg-V8-CT", 
  "NegCtrl-V16-CT"
)

reads_sum$description <- factor(reads_sum$description, levels = SAMPLE_ORDER)

# ── Log-transform reads ───────────────────────────────────
reads_sum <- reads_sum %>%
  mutate(reads_log = ifelse(reads_without_dups == 0, NA, log10(reads_without_dups)))

global_max <- max(reads_sum$reads_without_dups, na.rm = TRUE)
global_max_log <- log10(max(reads_sum$reads_without_dups, na.rm = TRUE))

# ── Heatmap Plot helper ───────────────────────────────────────────
make_heatmap <- function(df, title) {
  ggplot(df, aes(x = description, y = species, fill = reads_log)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "F", direction = -1, begin = 0.3, end = 0.8,
      na.value = "#bdd5e7",
      name = "Read counts\n(0 = light blue)",
      limits = c(0, global_max_log),
      breaks = c(0, 1, 2, 3, 4, 5, global_max_log),
      labels = c("1", "10", "100", "1000", "10000", "100000", as.character(round(global_max)))
    ) +
    geom_text(aes(label=reads_without_dups), color="black")+
    labs(title = title, x = "Sample", y = "Virus") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 10),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    )
}

heatmap_complete  <- make_heatmap(reads_sum,  "Virus detection across samples")

print(heatmap_complete)

ggsave(file.path(OUTPUT_DIR, "heatmap_complete.png"), heatmap_complete,  width = 14, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "heatmap_complete.svg"), heatmap_complete,  width = 14, height = 6, dpi = 300)

#---cDNA Synthesis impact-------------
CDNA_POS_NEG<-c(
  "Pos-V8-CT",
  "Pos-V20-CF",
  "Neg-V8-CT", 
  "Neg-V20-CF")

cdna_pos_neg_df<-reads_sum %>% 
  filter(description %in% CDNA_POS_NEG)

cdna_pos_neg_df$description<-factor(cdna_pos_neg_df$description, levels=CDNA_POS_NEG)

cDNA_alphatorque <- cdna_pos_neg_df %>%
  filter(species=="Alphatorquevirus")

cDNA_alphatorque_barplot <- ggplot(cDNA_alphatorque, aes(x=description, y=reads_without_dups))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1))+
  labs(title="Alphatorquvirus reads across samples\nAlignment, filtered, high coverage regions removed")

print(cDNA_alphatorque_barplot)
ggsave(file.path(OUTPUT_DIR,"cdna_pos_neg_anello_barplot.png"),  cDNA_alphatorque_barplot,  width = 6, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR,"cdna_pos_neg_anello_barplot.svg"),  cDNA_alphatorque_barplot,  width = 6, height = 6, dpi = 300)

#-------- Monkeypox Concentration Gradient -----
MONKEYPOX_CONC_ROW<-c(
  "Pos-V0-CF",
  "Pos-V2-CF", 
  "Pos-V20-CF",
  "Pos-V200-CF")

Monkeypox_conc_df<-reads_sum %>% 
  filter(description %in% MONKEYPOX_CONC_ROW)

Monkeypox_conc_df$description<-factor(Monkeypox_conc_df$description, levels=MONKEYPOX_CONC_ROW)

monkeypox_only_df <- Monkeypox_conc_df %>%
  filter(species=="Monkeypox_virus")

monkeypox_corr_plot <- ggplot(monkeypox_only_df, aes(x = mpox_gc_input, y = reads_without_dups)) +
  geom_point(shape = 4, size = 3, stroke = 1.3) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000, 10000)) +
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 2, 20, 200)) +
  labs(title = "Monkeypox reads vs Monkeypox gc input",
       x = "Monkeypox input (gc)", y = "Monkeypox reads") +
  theme_minimal()

print(monkeypox_corr_plot)
ggsave(file.path(OUTPUT_DIR, "monkeypox_corr.png"), monkeypox_corr_plot, width = 11, height = 6, dpi = 300) 
ggsave(file.path(OUTPUT_DIR, "monkeypox_corr.svg"), monkeypox_corr_plot, width = 11, height = 6, dpi = 300) 


#---EBV differences plot rpm-------------------------------------------------------------
EBV_only <- Monkeypox_conc_df %>%
  filter(species=="Human_gammaherpesvirus_4")

EBV_plot <- ggplot(EBV_only, aes(x=description, y=reads_without_dups))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1))+
  labs(title="EBV reads across samples\nAlignment, filtered, high coverage regions removed")

print(EBV_plot_rpm)
ggsave(file.path(OUTPUT_DIR, "EBV_barplot.png"), EBV_plot,  width = 6, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "EBV_barplot.svg"), EBV_plot,  width = 6, height = 6, dpi = 300)





