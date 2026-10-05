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
QPCR_REF <- "/path/to/Samples_overview_EBV_qPCR.xlsx" #file with additional sample information (e.g. Cp values)
VIRUS_LENGTH <- "/path/to/virus_length.txt" #for each virus its respective genome length, manually written

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
  "B01-Shearing-P4B",
  "B01-Shearing-P4C",
  "B01-Shearing-P4D",
  "B01-NoShearing-P2A",
  "B01-NoShearing-P2B",
  "B01-NoShearing-P2C",
  "B03-Shearing-P3D",
  "B03-Shearing-P3E",
  "B03-Shearing-P3F",
  "B03-NoShearing-P3A",
  "B03-NoShearing-P3B",
  "B03-NoShearing-P3C",
  "B04-Shearing-P4E",
  "B04-Shearing-P4F",
  "B04-Shearing-P4G",
  "B04-NoShearing-P2F",
  "B04-NoShearing-P2G",
  "B04-NoShearing-P2H",
  "B05-Shearing-P3G",
  "B05-Shearing-P3H",
  "B05-Shearing-P4A",
  "B05-NoShearing-P2D",
  "B05-NoShearing-P2E", 
  "Neg-NoShearing-P4H"
)

reads_sum$description <- factor(reads_sum$description, levels = SAMPLE_ORDER)
reads_sum$shearing <- !grepl(pattern = "NoShearing" , x = reads_sum$description)
reads_sum$sample_name <- stringr::str_sub(string=reads_sum$description, start=1, end=3)

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


#--viral load estimation------------------------------------------------

virus_length <- read_tsv(VIRUS_LENGTH)

reads_sum <- reads_sum %>%
  left_join(virus_length, by=c("species"="virus"))%>%
  group_by(description) %>%
  mutate(
    ref_reads = reads_without_dups[species == "Monkeypox_virus"],  
    reads_normalized = reads_without_dups / ref_reads, 
    estimated_gc=reads_normalized*196858*20/length,
    viral_load=estimated_gc*1000/3.8
  ) %>%
  ungroup() %>%
  select(-ref_reads)

reads_sum <- reads_sum %>%
  mutate(reads_normalized_log = ifelse(reads_normalized == 0, NA, log10(reads_normalized)), 
         estimated_gc_log = ifelse(estimated_gc == 0, NA, log10(estimated_gc)),
         viral_load_log = ifelse(viral_load == 0, NA, log10(viral_load)))


make_heatmap_viral_load <- function(df, title) {
  ggplot(df, aes(x = description, y = species, fill = viral_load_log))+
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "F", direction = -1, begin = 0.3, end = 0.8,
      na.value = "#bdd5e7",
      name = "viral load\n(no reads = light blue)",
      breaks = c(1, 2, 3, 4, 5, 6),
      labels = c("10", "100", "1k", "10k", "100k", "1M")
    ) +
    geom_text(aes(label=round(viral_load, 1)), color="black")+
    labs(title = title, x = "Sample", y = "Virus") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 10),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    )
}

heatmap_viral_load  <- make_heatmap_viral_load(reads_sum,  "estimated viral load in samples\nAlignment, filtered, high coverage regions removed")

print(heatmap_viral_load)

ggsave(file.path(OUTPUT_DIR, "heatmap_viral_load.png"), heatmap_viral_load,  width = 20, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "heatmap_viral_load.svg"), heatmap_viral_load,  width = 20, height = 6, dpi = 300)

#---EBV differences plots-------------------------------------------------------------
EBV_only <- reads_sum %>%
  filter(species=="Human_gammaherpesvirus_4")

EBV_plot_viral_load <- ggplot(EBV_only, aes(x=description, y=viral_load, fill=shearing))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.text.x  = element_text(size = 10, angle = 45, hjust = 1))+
  labs(title="viral load EBV across samples\nAlignment, filtered, high coverage regions removed")

print(EBV_plot_viral_load)
ggsave(file.path(OUTPUT_DIR, "EBV_barplot_viral_load.png"), EBV_plot_viral_load,  width = 14, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "EBV_barplot_viral_load.svg"), EBV_plot_viral_load,  width = 14, height = 6, dpi = 300)

#--- qPCR comparison ----------------------
qPCR_ref <- read_xlsx(QPCR_REF, sheet=3)
qPCR_comparison <- EBV_only %>%
  left_join(qPCR_ref, by=c("sample_name"="sample_number"))%>%
  group_by(sample_name)%>%
  mutate("mean_reads"=mean(reads_without_dups),  
         "mean_viral_load"=mean(viral_load),  
         "min_viral_load"  = min(viral_load),
         "max_viral_load"  = max(viral_load))

qPCR_comparison_plot_viral_load <- ggplot(qPCR_comparison)+
  geom_jitter(aes(x=Mean_Cp, y=viral_load, colour=shearing), size=5, alpha=0.5, height=50, width=0)+
  geom_point(aes(x=Mean_Cp, y=mean_viral_load), shape=4, size=6)+
  theme_minimal()

print(qPCR_comparison_plot_viral_load)
ggsave(file.path(OUTPUT_DIR, "qPCR_comparison_plot_viral_load.png"), qPCR_comparison_plot_viral_load,  width = 12, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "qPCR_comparison_plot_viral_load.svg"), qPCR_comparison_plot_viral_load,  width = 12, height = 6, dpi = 300)


#---- enrichment visualization---------------------------
qPCR_comparison_long <- qPCR_comparison %>%
  select(sample_name, mean_reads, EBV_read_counts_GS) %>%
  pivot_longer(cols = c(mean_reads, EBV_read_counts_GS),
               names_to = "method", values_to = "reads") %>%
  distinct()%>%
  mutate(method = recode(method,
                         enrichment_reads = "Targeted enrichment",
                         GS_reads = "Genome sequencing"))

enrichment_panel <- ggplot(qPCR_comparison_long, aes(x = sample_name, y = reads, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label=round(reads)), color="black",position = position_dodge(width = 0.7),
            vjust = -0.3)+
  scale_y_continuous(trans="pseudo_log", breaks=c(0,1,3,10,30,100,300)) +
  labs(y = "EBV reads (pseudo log scale)", x = "Sample", fill = NULL) +
  theme_minimal()

print(enrichment_panel)
ggsave(file.path(OUTPUT_DIR, "enrichment_panel.png"), enrichment_panel,  width = 6, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "enrichment_panel.svg"), enrichment_panel,  width = 6, height = 6, dpi = 300)

