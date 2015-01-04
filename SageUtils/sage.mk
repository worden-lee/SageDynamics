ifeq ($(HOME),)
export HOME = /usr/local/workingwiki/working-directories
export DOT_SAGE = $(HOME)/.sage
export SAGE_ROOT = /usr/share/sage
export LD_LIBRARY_PATH += :/usr/share/sage/local/lib
export SAGE_DOC = /usr/share/sage/devel/doc
endif

# how to run a sage script and capture its output
# the odd tee command is possibly temporary - it puts the 
# script's output into the .make.log (or console) as well as the output file
%.sage.out %.sage.tried : %.sage
	touch $*.sage.tried
	$(RM) $*.sage.out $*.sage.status $(STEP_PRODUCTS)
	(sage $< && touch $*.sage.status) | tee $*.sage.dmp
	[ -e $*.sage.status ] && ($(RM) $*.sage.status && mv $*.sage.dmp $*.sage.out) || exit 1

.PRECIOUS: %.sage.tried

# and a convenience rule to invoke sage interactively, but using all the
# exports from the projects' makefiles
sage :
	sage

sage-gdb :
	sage -gdb

.PHONY: sage sage-gdb
