DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: blast2matrix

blast2matrix: blast2matrix.00

blast2matrix.00:
	@../R_src/blast2matrix.r $@.blast8 $@.orfs2contigs.tsv > $@.out
	@$(CHECK)
