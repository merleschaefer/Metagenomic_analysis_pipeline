#--input and output - costumize these when using this script outside of pipeline------------------------------------------------------
input_coverage_data <- snakemake@input[["coverage"]]
input_fasta         <- snakemake@input[["fasta"]]
window_size         <- snakemake@params[["window_size"]]

out_csv   <- snakemake@output[["csv"]]
out_bed   <- snakemake@output[["bed"]]
out_txt   <- snakemake@output[["txt"]]
out_fasta <- snakemake@output[["fasta_out"]]
#--------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
#--------------------------------------------------------
#---Load coverage depth data-----------------------------
depth_df <- read_tsv(input_coverage_data, col_names = c("contig", "position", "depth"), col_types = cols(
  contig = col_character(),
  position = col_integer(),
  depth    = col_integer()
))
#---Sort into bins and flag------------------------------
bins_df <- depth_df |>
  mutate(bin = position %/% window_size) |>
  group_by(contig, bin) |>
  summarise(mean_depth = mean(depth), .groups = "drop") |>
  mutate(
    bin_start = bin * window_size + 1,
    bin_end   = bin_start + window_size - 1
  ) |>
  group_by(contig) |>
  mutate(
    q1        = quantile(mean_depth, 0.25),
    q3        = quantile(mean_depth, 0.75),
    iqr       = q3 - q1,
    threshold = 4 * iqr + q3 + 5,
    flagged   = mean_depth > threshold
  ) |>
  ungroup()

flagged_bins <- bins_df %>%
  filter(flagged) %>%
  arrange(contig, bin_start)

#---Fuse flagged bins into flagged regions, merge---------------
groups <- rep(1, nrow(flagged_bins))
if (nrow(flagged_bins) > 1) {
  for (i in 2:nrow(flagged_bins)) {
    same_virus <- flagged_bins$contig[i] == flagged_bins$contig[i - 1]
    adjacent   <- flagged_bins$bin_start[i] - 1 == flagged_bins$bin_end[i - 1]

    if (same_virus && adjacent) {
      groups[i] <- groups[i - 1]
    } else {
      groups[i] <- groups[i - 1] + 1
    }
  }
}
flagged_bins$group <- groups

merged <- flagged_bins %>%
  group_by(group, contig) %>%
  summarise(
    bin_start = min(bin_start),
    bin_end = max(bin_end),
    .groups = "drop"
  )

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write_csv(merged, out_csv)

#---write regions into BED file (0-based, half-open)------------------
bed_df <- merged %>%
  transmute(
    contig = contig,
    start  = bin_start - 1L,   # 1-based inclusive -> 0-based
    end    = bin_end           # numerically already correct as exclusive end
  ) %>%
  arrange(contig, start)
write_tsv(bed_df, out_bed, col_names = FALSE)

#---write regions into fasta file--------------------
regions <- sprintf("%s:%d-%d", merged$contig, merged$bin_start, merged$bin_end)
writeLines(regions, out_txt)

if (length(regions) > 0) {
  system2(
    "samtools",  # relies on samtools being on PATH (see environment.yml)
    c("faidx", input_fasta, "-r", out_txt),
    stdout = out_fasta
  )
} else {
  # No flagged regions at all: write an empty fasta so downstream rules
  # (which expect this file to exist) still succeed.
  file.create(out_fasta)
}
