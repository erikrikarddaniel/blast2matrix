#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

blast8file = '../test/blast2matrix.09.blast8'
orfs2contigsfile = '../test/blast2matrix.09.orfs2contigs.tsv'

blast = data.table(
    read_tsv(
    blast8file, 
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

orfs2contigs =  data.table(
  read_tsv(
    orfs2contigsfile,
    col_names = c('contig', 'seq'),
    col_types = c(contig = col_character(), query = col_character()),
    comment = '#'
  )
)

lengths = blast %>% filter(query==subject) %>%
  group_by(query) %>%
  summarise(len=sum(alignlen)) %>%
  ungroup() %>%
  select(query, len)

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
alldists = data.table(expand.grid(qcontig=unique(contigdists$qcontig), scontig=unique(contigdists$qcontig))) %>%
  left_join(contigdists, by=c('qcontig', 'scontig')) %>%
  replace_na(list('similarity'=0.0, 'one_minus_dist'=1.0))
