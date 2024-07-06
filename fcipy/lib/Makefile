# Build the ccq_py python module
MODULES := ci 

EXT_SUFFIX := $(shell python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
SOURCES := $(addsuffix .f90, $(MODULES))
TARGETS := $(addsuffix $(EXT_SUFFIX), $(MODULES))

# For GitHub CI pipeline, remove OpenMP
# the -mavx flag speeds things up!
FFLAGS := -O3 -ffree-line-length-512 -mavx -march=native -mtune=native# -fopenmp -lgomp

.PHONY: all clean

all: $(TARGETS)

%$(EXT_SUFFIX): %.f90
	echo $(MKLFOLDER)
	f2py --f90flags="$(FFLAGS)" -c $< -m $*

clean:
	rm $(TARGETS)
