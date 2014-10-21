ifeq ($(HOME),)
export HOME = /usr/local/workingwiki/working-directories
export DOT_SAGE = $(HOME)/.sage
export SAGE_ROOT = /usr/share/sage
export LD_LIBRARY_PATH += :/usr/share/sage/local/lib
export SAGE_DOC = /usr/share/sage/devel/doc
endif

# how to run a sage script and capture its output
# the odd tee command is possibly temporary - it puts the 
# script's output into the .make.log as well as the output file
%.sage.out : %.sage
	(sage $< | tee $<.dmp) && mv $<.dmp $@
