DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: blast2matrix

blast2matrix: blast2matrix.00 blast2matrix.01 blast2matrix.02 blast2matrix.03 blast2matrix.04 blast2matrix.05 blast2matrix.06

blast2matrix.00:
	@../R_src/blast2matrix --format=splitstree --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.01:
	@../R_src/blast2matrix --format=tsv --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.02:
	@../R_src/blast2matrix --symmetry=mean --format=tsv --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.03:
	@../R_src/blast2matrix --symmetry=max --format=tsv --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.04:
	@../R_src/blast2matrix --symmetry=min --format=tsv --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.05:
	@../R_src/blast2matrix --sim2dist=similarity --symmetry=mean --format=tsv --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)

blast2matrix.06:
	@../R_src/blast2matrix --sim2dist=similarity --symmetry=mean --format=abc --orfs2contigs=$@.orfs2contigs.tsv --blast8file=$@.blast8 > $@.out
	@$(CHECK)
