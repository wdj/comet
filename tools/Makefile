#------------------------------------------------------------------------------

sources := $(wildcard *.cc)

executables := $(patsubst %.cc,%,$(sources))

all_tools: $(executables)

test:	distclean all_tools
	./test.sh

clean distclean:
	rm -f $(executables)

.PHONY: all clean distclean
#------------------------------------------------------------------------------
