#!/usr/bin/env Rscript

# blast2matrix.r
#
# Reads a file with the "-m 8" tabular BLAST output format (LAST's version is
# also fine) and outputs a distance matrix, either from the entries in the
# input file or from some mapping of those to something else (e.g. ORF to
# contig).
#
# Author: daniel.lundin@lnu.se

# Required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

FORMAT_SPLITSTREE = 'splitstree'
FORMAT_TSV        = 'tsv'
FORMAT_DEFAULT    = FORMAT_TSV
FORMATS           = c(FORMAT_SPLITSTREE, FORMAT_TSV)

# Get arguments
option_list = list(
  make_option(
    c('--blast8file'), type='character',
    help='Name of file BLAST tabular file (-m 8); LAST\'s output also works'
  ),
  make_option(
    c('--format'), type='character', default=FORMAT_DEFAULT,
    help='Output format, default: [default]'
  ),
  make_option(
    c('--formats'), action='store_true', default=FALSE,
    help='Lists supported formats'
  ),
  make_option(
    c('--orfs2contigs'), type='character',
    help='Name of file with contig to orf mapping; format: contig<tab>orf'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$formats ) {
  write(cat("Supported formats: ", FORMATS, "\n"))
  quit('no')
}

logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

logmsg(sprintf("Reading blast8 file %s", opt$blast8file))
blast = data.table(
    read_tsv(
    opt$blast8file, 
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
logmsg(sprintf("Read %d sequences, %d hits", length(unique(blast$query)), length(blast$query)))

logmsg(sprintf("Reading opt$orfs2contigs file %s", opt$orfs2contigs)) 
o2f =  data.table(
  read_tsv(
    opt$orfs2contigs,
    col_names = c('contig', 'seq'),
    col_types = c(contig = col_character(), query = col_character()),
    comment = '#'
  )
)

logmsg("Calculating lengths")
lengths = blast %>% filter(query==subject) %>%
  group_by(query) %>%
  summarise(len=sum(alignlen)) %>%
  ungroup() %>%
  select(query, len)

logmsg("Calculating distances between *sequences*")
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

logmsg("Calculating distances between *contigs*")
contigdists = data.table(
  distances %>%
    inner_join(o2f %>% select(qcontig = contig, query = seq), by = c('query')) %>%
    inner_join(o2f %>% select(scontig = contig, subject = seq), by = c('subject')) %>%
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

logmsg(sprintf("Writing result in % s format", opt$format))

if ( opt$format == FORMAT_SPLITSTREE ) {
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
} else {
  logmsg(sprintf("Format '%s' not supported, see --formats for supported formats", opt$format, llevel = 'ERROR'))
  quit('no')
}

logmsg("Done")