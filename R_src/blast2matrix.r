#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(trailingOnly = T)

write(sprintf("%s: Reading blast8 file %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), args[1]), stderr())
blast = data.table(
    read_tsv(
    args[1], 
    col_names=c(
      'query', 'subject', 'percid', 'alignlen', 'mismatches', 'gapopen', 
      'qstart', 'qend', 'sstart', 'send', 'e_value', 'bitscore'
    ),
    col_types = c(
      query = col_character(), subject = col_character(),
      percid = col_double(), alignlen = col_integer(), mismatches = col_integer(),
      gapopen = col_integer(), 
      qstart = col_integer(), qend = col_integer(),
      sstart = col_integer(), send = col_integer(),
      e_value = col_double(), bitscore = col_double()
    ),
    comment = '#'
  ) %>%
    transmute(query, subject, percid, alignlen, n_id=percid/100.0*alignlen),
  key = c('query', 'subject')
)
write(sprintf("%s: Read %d sequences, %d hits", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(unique(blast$query)), length(blast$query)), stderr())

write(sprintf("%s: Reading orfs2contigs file %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), args[2]), stderr())
orfs2contigs =  data.table(
  read_tsv(
    args[2],
    col_names = c('contig', 'seq'),
    col_types = c(contig = col_character(), query = col_character()),
    comment = '#'
  )
)

write(sprintf("%s: Calculating lengths", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stderr())
lengths = blast %>% filter(query==subject) %>%
  group_by(query) %>%
  summarise(len=sum(alignlen)) %>%
  ungroup() %>%
  select(query, len)

write(sprintf("%s: Calculating distances between sequences", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stderr())
distances = data.table(
  blast %>%
    group_by(query, subject) %>%
    summarise(n_id = sum(n_id)) %>%
    ungroup() %>%
    inner_join(lengths, by='query') %>%
    transmute(query, subject, similarity = n_id/len) %>%
    mutate(one_minus_dist = ( 1 - similarity )),
  key = c('query', 'subject')
)

write(sprintf("%s: Calculating distances between contigs", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stderr())
contigdists = data.table(
  distances %>%
    inner_join(orfs2contigs %>% select(qcontig = contig, query = seq), by = c('query')) %>%
    inner_join(orfs2contigs %>% select(scontig = contig, subject = seq), by = c('subject')) %>%
    group_by(qcontig, scontig, query) %>%
    summarise(similarity=max(similarity), one_minus_dist=min(one_minus_dist)) %>%
    ungroup() %>%
    group_by(qcontig, scontig) %>%
    summarise(
      similarity = mean(similarity),
      one_minus_dist = mean(one_minus_dist)
    ) %>%
    ungroup(),
  key = c('qcontig', 'scontig')
)

# Make sure we get distances between all pairs
contigdists = data.table(expand.grid(qcontig=unique(contigdists$qcontig), scontig=unique(contigdists$qcontig))) %>%
  left_join(contigdists, by=c('qcontig', 'scontig')) %>%
  replace_na(list('similarity'=0.0, 'one_minus_dist'=1.0))

widedists = contigdists %>% select(-similarity) %>%
  spread(scontig, one_minus_dist, fill=1.0)

write(sprintf("%s: Writing result to file", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stderr())
cat("#nexus\n")
cat("BEGIN Taxa;\n")
cat(sprintf("DIMENSIONS ntax=%d;\n", length(widedists$qcontig)))
cat("TAXLABELS\n")
cat(sprintf("[%d] '%s'", 1:length(widedists$qcontig), widedists$qcontig), sep="\n")
cat(";\n")
cat("END; [Taxa]\n")
cat("BEGIN Distances;\n")
cat(sprintf("DIMENSIONS ntax=%d;\n", length(widedists$qcontig)))
cat("FORMAT labels=no diagonal triangle=both;\n")
cat("MATRIX\n")
write.table(widedists %>% select(-qcontig), sep="\t", row.names=F, col.names=F)
cat(";\n")
cat("END; [Distances]\n")

write(sprintf("%s: Done", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stderr())
