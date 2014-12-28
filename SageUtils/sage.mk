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
%.sage.out : %.sage
	$(RM) $@ $*.status
	(sage $< && touch $*.status) | tee $*.dmp
	[ -e $*.status ] && ($(RM) $*.status && mv $*.dmp $@) || exit 1

#	(sage $< | tee $<.dmp) && mv $<.dmp $@

# and a convenience rule to invoke sage interactively, but using all the
# exports from the projects' makefiles
sage :
	sage

sage-gdb :
	sage -gdb

.PHONY: sage
