# Compile the custom Grid measurement programs against the installed libGrid.
# Uses grid-config (installed by the build under $PREFIX/bin) for the exact nvcc flags the
# Grid build used -- so our programs match the library's SIMD/accelerator configuration.
#
#   make GRID_CONFIG=$HOME/ed-lattice/bin/grid-config
PREFIX      ?= $(HOME)/ed-lattice
GRID_CONFIG ?= $(PREFIX)/bin/grid-config

CXX      := $(shell $(GRID_CONFIG) --cxx 2>/dev/null)
CXXFLAGS := $(shell $(GRID_CONFIG) --cxxflags 2>/dev/null)
LDFLAGS  := $(shell $(GRID_CONFIG) --ldflags 2>/dev/null)
LIBS     := $(shell $(GRID_CONFIG) --libs 2>/dev/null)

PROGS := measure_potential generate_gauge measure_baryon measure_baryon_gevp generate_dynamical measure_modenumber measure_bag_profile generate_sextet

all: check $(PROGS)

check:
	@if [ -z "$(CXX)" ]; then \
	  echo "ERROR: grid-config not found at $(GRID_CONFIG)."; \
	  echo "Set GRID_CONFIG=<grid prefix>/bin/grid-config (default \$$PREFIX/bin)."; \
	  exit 1; fi
	@echo "using CXX=$(CXX)"

%: %.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	rm -f $(PROGS)

.PHONY: all check clean
