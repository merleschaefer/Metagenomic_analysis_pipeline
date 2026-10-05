library(dplyr)
library(readr)
library(writexl)

insert_sizes <- read.table("/path/to/insert_sizes.txt", header=TRUE) #insert sizes of final reads
USI_description <- read.csv("/path/to/USI_description.csv") #for translation of lab internal codes to sample descriptions for interpretability
insert_sizes_anello <- read.table("/path/to/insert_sizes_anello.txt", header=TRUE) #insert sizes of final anellovirus reads
total_reads <- read.csv("/path/to/total_read_counts.csv") #total number of reads per sample before host removal
reads_without_host <- read.csv("/path/to/host_removed_read_counts.csv")  #total number of reads per sample after host removal
read_counts <- read_tsv("/path/to/final_read_counts.tsv") #final read counts with, without duplicates

output <- "path/to/test_run_2_stats.xlsx"

#insert sizes
insert_sizes <- insert_sizes %>%
  filter(sample!="merged_filtered")

insert_sizes$sample <- sub("_.*", "", insert_sizes$sample)

insert_sizes <- insert_sizes%>%
  left_join(USI_description, by=c("sample"="USI"))

insert_sizes <- insert_sizes %>%
  filter(size != "NA") %>%
  group_by(sample) %>%
  mutate(mean_insert_size = mean(size))

insert_sizes$shearing <- !grepl(pattern = "*NoShearing*" , x = insert_sizes$description)
insert_sizes$sample_name <- stringr::str_sub(string=insert_sizes$description, start=1, end=3)

insert_sizes_means <-insert_sizes%>%
  select(c(-size)) %>%
  unique()

#insert sizes anelloviruses

insert_sizes_anello <- insert_sizes_anello %>%
  filter(sample!="merged_filtered")

insert_sizes_anello$sample <- sub("_.*", "", insert_sizes_anello$sample)

insert_sizes_anello <- insert_sizes_anello%>%
  left_join(USI_description, by=c("sample"="USI"))

insert_sizes_anello <- insert_sizes_anello %>%
  filter(size != "NA") %>%
  group_by(sample) %>%
  mutate(mean_insert_size_anello = mean(size))

insert_sizes_anello_means <-insert_sizes_anello%>%
  select(c(-size)) %>%
  unique()

# host percentage
host_percentages <- total_reads %>%
  inner_join(reads_without_host, by="sample", suffix = c("_with_host", "_without_host"))
host_percentages$sample <- sub("_.*", "", host_percentages$sample)
host_percentages <- host_percentages%>%
  inner_join(USI_description, by=c("sample"="USI"))

host_percentages$host_percentage <- (((host_percentages$r1_reads_with_host)-(host_percentages$r1_reads_without_host))*100)/(host_percentages$r1_reads_with_host)

#duplicates
read_counts <- read_counts %>%
  filter(sample!="merged_filtered.clean")

read_counts$sample <- sub("_.*", "", read_counts$sample)
duplicates <-read_counts%>%
  filter(reads_with_dups>0)%>%
  group_by(sample)%>%
  mutate(percent_duplicates=(sum(reads_with_dups)-sum(reads_without_dups))*100/sum(reads_with_dups))%>%
  select(sample, percent_duplicates)%>%
  unique()

#combine everything
df <- host_percentages%>%
  left_join(insert_sizes_means, by=c("sample","description"))%>%
  left_join(insert_sizes_anello_means, by=c("sample","description"))%>%
  left_join(duplicates, by="sample")%>%
  select(c("description","host_percentage","mean_insert_size", "mean_insert_size_anello", "percent_duplicates"))

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

df$description <- factor(df$description, levels = SAMPLE_ORDER)
df <- df %>% arrange(description)

write_xlsx(df, output)
