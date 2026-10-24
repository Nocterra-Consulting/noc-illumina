# Makefile for the ILLUMINA Fortran kernel.
#
# Targets:
#   all     (default) release build -> bin/illumina, bin/continue_illumina
#   debug   checked build           -> bin/illumina_debug
#   openmp  release build + OpenMP  -> bin/illumina_omp
#   pgo     release build with link-time and profile-guided optimisation
#           -> bin/illumina, bin/continue_illumina; the profile comes
#           from a short synthetic run (about 4 s). The results are the
#           same as the release build (checked with gfortran 13 and 14);
#           run "make test" after it on a new compiler.
#   test    build "all", then run the synthetic regression harness
#           (case_small, case_hill, case_cloud_hill, case_cloud_2nd)
#           (also runs test-omp when bin/illumina_omp exists)
#   test-omp build "openmp", then check it against the serial build
#   clean   remove the build directory and the binaries
#
# The kernel source is fixed-form Fortran 77 and is not preprocessed.
# The version string is therefore substituted with sed in a copy of
# illumina.f inside $(BUILD). The files in kernel/ are never modified.

# make pre-defines FC=f77; override that default but respect a user setting.
ifeq ($(origin FC),default)
FC = gfortran
endif
FFLAGS   ?= -O3 -mcmodel=medium -Wunused-parameter
DBGFLAGS ?= -O0 -g -mcmodel=medium -Wall -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow
# -fopenmp implies -frecursive, which moves every local array to the
# stack. zondif(3000000,3) in the main program is 36 MB and overflows the
# default 8 MB stack. The explicit size limit keeps the default placement:
# arrays above the limit stay static (main program only), everything
# else is on the stack and thread-safe.
OMPFLAGS ?= $(FFLAGS) -fopenmp -fmax-stack-var-size=65536
PYTHON   ?= python3

KERNEL := kernel
BUILD  := build
BIN    := bin

VERSION := $(shell sed -n 's/^__version__ *= *"\(.*\)"/\1/p' illum/__init__.py)

# Explicit list (mirrors bin/makeILLUMINA). continue_illumina.f, ies2fctem.f,
# cloudtransmitance.f and FortranTools/ are intentionally excluded.
KERNEL_SRCS := \
	zone_diffusion.f \
	diffusion.f \
	angle3points.f \
	anglesolide.f \
	transmita.f \
	transmitm.f \
	transmitl.f \
	anglezenithal.f \
	angleazimutal.f \
	horizon.f \
	cloudreflectance.f \
	twodin.f \
	twodout.f \
	transTOA.f \
	curvature.f \
	molabs.f \
	expzh.f

MAIN_SRC := $(BUILD)/illumina.f
SRCS     := $(MAIN_SRC) $(addprefix $(KERNEL)/,$(KERNEL_SRCS))

.PHONY: all debug openmp pgo test test-omp clean

all: $(BIN)/illumina $(BIN)/continue_illumina

debug: $(BIN)/illumina_debug

openmp: $(BIN)/illumina_omp

$(BUILD) $(BIN):
	mkdir -p $@

# Copy of the main program with the version string substituted.
$(MAIN_SRC): $(KERNEL)/illumina.f illum/__init__.py | $(BUILD)
	sed 's/__version__/$(VERSION)/' $< > $@

$(BIN)/illumina: $(SRCS) | $(BIN)
	$(FC) $(FFLAGS) $(SRCS) -o $@

$(BIN)/illumina_debug: $(SRCS) | $(BIN)
	$(FC) $(DBGFLAGS) $(SRCS) -o $@

$(BIN)/illumina_omp: $(SRCS) | $(BIN)
	$(FC) $(OMPFLAGS) $(SRCS) -o $@

# The profile file names follow the output path, so both passes write
# $(BIN)/illumina.
PGO_DIR  := $(abspath $(BUILD)/pgo)
PGO_CASE := $(BUILD)/pgo_case
pgo: $(SRCS) $(BIN)/continue_illumina | $(BIN)
	rm -rf $(PGO_DIR) $(PGO_CASE)
	$(FC) $(FFLAGS) -flto=auto -fprofile-generate=$(PGO_DIR) $(SRCS) -o $(BIN)/illumina
	$(PYTHON) tests/regression/make_case.py --out $(PGO_CASE) --size 101 --dx 50 --hill 150 600 1200 800 --hill-lamps > /dev/null
	printf '10 30\n30 120\n60 250\n5 300\n' > $(PGO_CASE)/angles.txt
	cd $(PGO_CASE) && $(abspath $(BIN))/illumina illumina.in synth angles.txt nomaps nolog > /dev/null
	$(FC) $(FFLAGS) -flto=auto -fprofile-use=$(PGO_DIR) -fprofile-partial-training $(SRCS) -o $(BIN)/illumina

$(BIN)/continue_illumina: $(KERNEL)/continue_illumina.f | $(BIN)
	$(FC) $< -o $@

test: all
	$(PYTHON) tests/regression/run_regression.py --binary $(BIN)/illumina --case tests/regression/case_small
	$(PYTHON) tests/regression/run_regression.py --binary $(BIN)/illumina --case tests/regression/case_hill
	$(PYTHON) tests/regression/run_regression.py --binary $(BIN)/illumina --case tests/regression/case_cloud_hill
	$(PYTHON) tests/regression/run_regression.py --binary $(BIN)/illumina --case tests/regression/case_cloud_2nd
	$(PYTHON) tests/regression/run_angles.py --binary $(BIN)/illumina --case tests/regression/case_small
	@if [ -x $(BIN)/illumina_omp ]; then $(MAKE) --no-print-directory test-omp; fi

# OpenMP check: 1 thread must be bit-identical to the serial build,
# 4 threads must agree within rtol 1e-4.
test-omp: all openmp
	$(PYTHON) tests/regression/run_omp.py --serial $(BIN)/illumina --omp $(BIN)/illumina_omp --case tests/regression/case_small
	$(PYTHON) tests/regression/run_omp.py --serial $(BIN)/illumina --omp $(BIN)/illumina_omp --case tests/regression/case_hill
	$(PYTHON) tests/regression/run_omp.py --serial $(BIN)/illumina --omp $(BIN)/illumina_omp --case tests/regression/case_cloud_hill
	$(PYTHON) tests/regression/run_omp.py --serial $(BIN)/illumina --omp $(BIN)/illumina_omp --case tests/regression/case_cloud_2nd

clean:
	rm -rf $(BUILD)
	rm -f $(BIN)/illumina $(BIN)/illumina_debug $(BIN)/illumina_omp $(BIN)/continue_illumina
