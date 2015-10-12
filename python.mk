%.py.out : %.py
	(python $< >$<.dmp && mv $<.dmp $@) || ! cat $<.dmp
