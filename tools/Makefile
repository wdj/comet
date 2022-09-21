#------------------------------------------------------------------------------

sources := $(wildcard *.cc)

executables := $(patsubst %.cc,%,$(sources))

all_tools: $(executables)

clean distclean:
	rm -f $(executables)

.PHONY: all clean distclean
#------------------------------------------------------------------------------
