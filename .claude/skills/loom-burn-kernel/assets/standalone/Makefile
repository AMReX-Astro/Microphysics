# Standalone burn reproducer (see README.md). Needs only ROCm's hipcc.
#
#   make                      # MATH=ocml: hipcc's own device math (as the stock HIP build)
#   make MATH=recipes         # the Loom path's f64 recipes (f64_math*.H) instead
#   make ARCH=gfx942          # GPU target (default: first GPU amdgpu-arch reports, else gfx90a)
#   make asm                  # device assembly of the kernel (burn_standalone_$(MATH).s)
#
# -ffp-contract=off is required for bit-exact comparison with the stock
# Microphysics reference build (which uses it too); keep it.

HIPCC ?= hipcc
ARCH  ?= $(shell amdgpu-arch 2>/dev/null | head -1)
ifeq ($(ARCH),)
  ARCH := gfx90a
endif
MATH  ?= ocml
EXTRA ?=

CXXFLAGS := -std=c++20 -O3 -ffp-contract=off --offload-arch=$(ARCH) -DMPL_HIP_STANDALONE -I.
ifeq ($(MATH),recipes)
  CXXFLAGS += -DMPL_HIP_MATH_RECIPES -Ishim
else ifneq ($(MATH),ocml)
  $(error MATH must be ocml or recipes)
endif

EXE     := burn_standalone_$(MATH)
SOURCES := main.hip burn_kernel.cpp $(wildcard *.H) $(wildcard shim/loomcxx/*.h)

$(EXE): $(SOURCES)
	$(HIPCC) $(CXXFLAGS) $(EXTRA) main.hip -o $@

asm: $(EXE).s
$(EXE).s: $(SOURCES)
	$(HIPCC) $(CXXFLAGS) $(EXTRA) --offload-device-only -S main.hip -o $@ 2>/dev/null || \
	$(HIPCC) $(CXXFLAGS) $(EXTRA) --cuda-device-only -S main.hip -o $@

clean:
	rm -f burn_standalone_ocml burn_standalone_recipes burn_standalone_*.s

.PHONY: asm clean
