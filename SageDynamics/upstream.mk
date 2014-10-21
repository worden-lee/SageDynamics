# The basic makefile stuff is upstream, import it from that project to this one
-include $(SageUtils)/sage.mk
-include $(SageUtils)/step.mk

# When we need something from an upstream project, make it there 
$(SageUtils)/% :
	$(MAKE) -C $(SageUtils) $*
