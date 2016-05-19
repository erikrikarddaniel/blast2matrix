DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: blast2matrix

blast2matrix: blast2matrix.00 blast2matrix.01 blast2matrix.02 blast2matrix.03 blast2matrix.04 blast2matrix.05 blast2matrix.06 blast2matrix.07

blast2matrix.00:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=tsv $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.01:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=tsv --symmetry=mean $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.02:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=tsv --symmetry=max $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.03:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=tsv --symmetry=min $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.04:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=phylip --symmetry=min $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.05:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=splitstree --symmetry=min $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.06:
	@../ruby_src/blast2matrix --distance=orfs --sim2dist=one_minus --format=tsv --symmetry=mean $@.blast8 > $@.out
	@$(CHECK)

blast2matrix.07:
	@../ruby_src/blast2matrix --distance=contigs --orfs2contigs=$@.orfs2contigs.tsv --sim2dist=one_minus --format=tsv --symmetry=mean $@.blast8 > $@.out
	@$(CHECK)
