# Makefile to connect the Sage code here with the code in SageDynamics
 
# When we need something from an upstream project, make it there
$(SageDynamics)/% :
	$(MAKE) -C $(SageDynamics) $*

$(SageUtils)/% :
	$(MAKE) -C $(SageUtils) $*
 
# The good makefile stuff is upstream, just reuse it
-include $(SageUtils)/sage.mk
-include $(SageUtils)/step.mk
