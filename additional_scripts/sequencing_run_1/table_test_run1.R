library(dplyr)
library(readr)
library(writexl)

insert_sizes <- read.table("/path/to/insert_sizes.txt", header=TRUE) #insert sizes of final reads
USI_description <- read.csv("/path/to/USI_description.csv") #for translation of lab internal codes to sample descriptions for interpretability
total_reads <- read.csv("/path/to/total_read_counts.csv") #total number of reads per sample before host removal
reads_without_host <- read.csv("/path/to/host_removed_read_counts.csv")  #total number of reads per sample after host removal
read_counts <- read_tsv("/path/to/final_read_counts.tsv") #final read counts with, without duplicates

output <- "/path/to/output/test_run_1_stats.xlsx"

#insert size
insert_sizes <- insert_sizes %>%
  filter(sample!="merged_filtered")

insert_sizes$sample <- sub("_.*", "", insert_sizes$sample)

insert_sizes <- insert_sizes%>%
  left_join(USI_description, by=c("sample"="USI"))

insert_sizes <- insert_sizes %>%
  filter(size != "NA") %>%
  group_by(sample) %>%
  mutate(mean_insert_size = mean(size))


insert_sizes_means <-insert_sizes%>%
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
  left_join(duplicates, by="sample")%>%
  select(c("description","host_percentage","mean_insert_size","percent_duplicates"))
write_xlsx(df, output)

