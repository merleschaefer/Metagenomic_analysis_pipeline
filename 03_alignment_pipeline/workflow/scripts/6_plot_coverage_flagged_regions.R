#--input and output - costumize these when using this script outside of pipeline------------------------------------------------------
input_coverage_data  <- snakemake@input[["coverage"]]
input_contig_names   <- snakemake@input[["contig_names"]]
input_flagged_regions <- snakemake@input[["flagged_csv"]]

out_png <- snakemake@output[["png"]]
out_svg <- snakemake@output[["svg"]]
#--------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
#---Load data-----------------------------
depth_df <- read_tsv(input_coverage_data, col_names = c("contig", "position", "depth"), col_types = cols(
  contig = col_character(),
  position = col_integer(),
  depth    = col_integer()
))
flagged_150 <- read.csv(input_flagged_regions)
contig_names <- read.csv(input_contig_names)

#---merge with name of virus------------------
depth_df <- depth_df %>%
  inner_join(contig_names, by = c("contig" = "identifier"))
flagged_150 <- flagged_150 %>%
  inner_join(contig_names, by = c("contig" = "identifier"))

#filter for min 2 depth
filtered_df <- depth_df %>%
  group_by(contig) %>%
  filter(max(depth) >= 2) %>%
  ungroup()

#calculate max depth per contig
max_depth_df <- depth_df %>%
  group_by(contig, species) %>%
  summarise(max_depth = max(depth), .groups = "drop")

#make bottom tracks for flaggings
add_bottom_track <- function(df, track_num, total_tracks = 4) {
  df %>%
    left_join(max_depth_df, by = c("contig", "species")) %>%
    mutate(
      track_height = max_depth * 0.08,   # each track = 8% of max depth
      gap = max_depth * 0.02,            # spacing between tracks
      ymax = -(track_num - 1) * (track_height + gap),
      ymin = ymax - track_height
    )
}
flagged_150 <- add_bottom_track(flagged_150, 1)

#---plot----------------------------------------
plot <- ggplot(filtered_df) +
  geom_line(aes(x = position, y = depth)) +
  geom_rect(
    data = flagged_150,
    aes(xmin = bin_start, xmax = bin_end,
        ymin = ymin, ymax = ymax,
        fill = "Window size 150"),
    inherit.aes = FALSE
  ) +
  facet_wrap(~species, scales = "free") +
  scale_fill_manual(
    name = "Flagged bins",
    values = c("Window size 150" = "red", "Window size 500" = "blue")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_minimal()

ggsave(out_png, plot, width = 25, height = 10, dpi = 100)
ggsave(out_svg, plot, width = 25, height = 10, dpi = 100)
