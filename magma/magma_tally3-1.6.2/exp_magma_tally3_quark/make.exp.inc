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

LIB += -L$(MAGMA_tally3_DIR)/lib -L$(MAGMA_tally3_DIR)/exp_magma_tally3_quark -lmagma_tally3_quark

#
# Quark
#
QUARKDIR ?= $(MAGMA_tally3_DIR)/quark
INC += -I$(QUARKDIR)
LIB += -L$(QUARKDIR) -lquark

INC += $(MORE_INC) -I../include -I../schedule -I../control -I../core -I../debug -I../magma_tally3_insert -I../magma_tally3_task

LIBMAGMA_tally3_QUARK = $(MAGMA_tally3_DIR)/exp_magma_tally3_quark/libmagma_tally3_quark.a

#
# 
#
INC += -I../../testing
LIB += -L../../testing -ltest
