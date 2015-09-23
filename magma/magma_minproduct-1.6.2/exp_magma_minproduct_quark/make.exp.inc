MORE_INC=
#
# CPU peak performance (for one cpu only)
# 
# Leave blank is unknown
# CPU_PEAK = 11.2
#MORE_INC += -DCPU_PEAK=8.4 

#
# GPU peak performance (for one GPU only)
# 
# Leave blank is unknown
# GPU_PEAK = 504.0
#MORE_INC += -DGPU_PEAK=504.0

LIB += -L$(MAGMA_minproduct_DIR)/lib -L$(MAGMA_minproduct_DIR)/exp_magma_minproduct_quark -lmagma_minproduct_quark

#
# Quark
#
QUARKDIR ?= $(MAGMA_minproduct_DIR)/quark
INC += -I$(QUARKDIR)
LIB += -L$(QUARKDIR) -lquark

INC += $(MORE_INC) -I../include -I../schedule -I../control -I../core -I../debug -I../magma_minproduct_insert -I../magma_minproduct_task

LIBMAGMA_minproduct_QUARK = $(MAGMA_minproduct_DIR)/exp_magma_minproduct_quark/libmagma_minproduct_quark.a

#
# 
#
INC += -I../../testing
LIB += -L../../testing -ltest
