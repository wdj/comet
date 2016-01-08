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

LIB += -L$(MAGMA_tally4_DIR)/lib -L$(MAGMA_tally4_DIR)/exp_magma_tally4_quark -lmagma_tally4_quark

#
# Quark
#
QUARKDIR ?= $(MAGMA_tally4_DIR)/quark
INC += -I$(QUARKDIR)
LIB += -L$(QUARKDIR) -lquark

INC += $(MORE_INC) -I../include -I../schedule -I../control -I../core -I../debug -I../magma_tally4_insert -I../magma_tally4_task

LIBMAGMA_tally4_QUARK = $(MAGMA_tally4_DIR)/exp_magma_tally4_quark/libmagma_tally4_quark.a

#
# 
#
INC += -I../../testing
LIB += -L../../testing -ltest
