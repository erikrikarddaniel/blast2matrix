# blast2matrix

blast2matrix converts a blast8 formated input file into a pairwise distance file.

The program is an R script located in the `R_src` directory (there's also an old Ruby
implementation (hopelessly slow on large data sets) in the `ruby_src` subdirectory).

Distances are based on fraction identity between pairs, seen over the length of the query.
Since fraction identity is a similarity score, distances are created by taking one minus
the similarity. The script is open for other implementations.

Not implemented in the R script yet, but in the Ruby version: The program can convert 
distances to make the matrix symmetrical by taking the mean, minimum or maximum of 
distances between pairs.

The R version can only output distances between contigs, while the Ruby version can also
convert distances between individual entities (e.g. proteins or orfs) to distances between
higher level entities (e.g. genomes or contigs). Currently this is calculated as the mean 
of distances less than a threshold (1.0 by default). The threshold can be set on the command
line.

Currently there are two output formats (three in the Ruby implementation): tab separated
(, phylip) and splitstree (Nexus).

See `--help` for instructions.
