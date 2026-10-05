# ── Libraries ────────────────────────────────────────────
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
# ── Paths ────────────────────────────────────────────────
BRACKEN_FILE <- "/path/to/pipeline_results/results/bracken/bracken_minusb_combined_reports.txt"
VIRAL_TAXIDS <- "./viral_taxids.txt"
USI_description <- "/path/to/USI_description.csv" # a file connecting the USI from sequencing to a description of the sample to make analyses easier

OUTPUT_DIR <- "/path/to/output"


# ── Load bracken results ────────────────────────────────────
bracken<- read.table(file=BRACKEN_FILE, sep = '\t',header = TRUE)

# ── Keep only _num columns + identifiers ─────────────────
num_cols <- c("name", "taxonomy_id", colnames(bracken)[grepl("_num$", colnames(bracken))])
bracken <- bracken %>%
  select(all_of(num_cols))

#  Trim DE column names at the first underscore 
colnames(bracken) <- ifelse(
  startsWith(colnames(bracken), "DE"),
  sub("_.*", "", colnames(bracken)),
  colnames(bracken)
)

# Replace DE column names with sample descriptions from CSV
lookup <- read.csv(USI_description, stringsAsFactors = FALSE)

# Build a named vector: names = USI, values = description
name_map <- setNames(lookup$description, lookup$USI)

# Replace column names where a match exists, leave others unchanged
colnames(bracken) <- ifelse(
  colnames(bracken) %in% names(name_map),
  name_map[colnames(bracken)],
  colnames(bracken)
)


# read the list of viral taxids (one taxid per line, no header)
viral_taxids <- read_lines(VIRAL_TAXIDS) %>% as.integer()

# filter to keep only rows whose taxid is in the viral list
bracken_viruses <- bracken %>%
  filter(taxonomy_id %in% viral_taxids)

# pivot to long format
bracken_viruses_long <- bracken_viruses %>%
  pivot_longer(cols = c(-name, -taxonomy_id), names_to = "description", values_to = "reads") 

#order viruses by prevalence
order_viruses <- function(df) {
  virus_order <- df %>%
    group_by(name) %>%
    summarise(n_detected = sum(reads > 0)) %>%
    arrange(desc(n_detected)) %>%
    pull(name)
  
  virus_levels <- c(rev(virus_order))
  
  df <- df %>%
    mutate(
      name = factor(name, levels = virus_levels)
    )
  }

bracken_viruses_long  <- order_viruses(bracken_viruses_long)


# ── Log-transform reads ───────────────────────────────────
bracken_viruses_long <- bracken_viruses_long %>%
  mutate(reads_log = ifelse(reads == 0, NA, log10(reads)))

#adjust sample order
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

bracken_viruses_long$description <- factor(bracken_viruses_long$description, levels = SAMPLE_ORDER)

# ── Heatmap Plot helper ───────────────────────────────────────────

make_heatmap <- function(df, title) {
  ggplot(df, aes(x = description, y = name, fill = reads_log)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "F", direction = -1, begin = 0.3, end = 0.8,
      na.value = "#bdd5e7",
      name = "Read counts\n(0 = light blue)",
      limits = c(0, 7),
      breaks = c(0, 1, 2, 3, 4, 5,6,7),
      labels = c("1", "10", "100", "1k", "10k", "100k", "1M", "10M")
    ) +
    geom_text(aes(label=reads), color="black")+
    labs(title = title, x = "Sample", y = "Virus") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 10),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    )
}

heatmap_complete <- make_heatmap(bracken_viruses_long,  "Virus detection across samples")

print(heatmap_complete)

ggsave(file.path(OUTPUT_DIR, "bracken_heatmap_complete.png"), heatmap_complete,  width = 10, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "bracken_heatmap_complete.svg"), heatmap_complete,  width = 10, height = 6, dpi = 300)



