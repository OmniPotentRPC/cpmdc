# Plane-wave scaffold. It does not link the embed library.
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
MPICC ?= mpicc
CFLAGS := -O2 -std=c11 -Wall -Wextra -Wpedantic -ffp-contract=off -fopenmp
CPPFLAGS := -I$(ROOT)/include -I$(ROOT)/bench -DCPMDC_PWSCALE_HAVE_MPI=1 -DCPMDC_PWSCALE_HAVE_FFTW=1
LDFLAGS := -fopenmp
LIBS := -lm -lfftw3_omp -lfftw3

BUILD := $(ROOT)/bench/build

.PHONY: all test clean

all: $(BUILD)/cpmdc-pwscale $(BUILD)/test_pwscale_agreement \
     $(BUILD)/test_pwscale_collectives $(BUILD)/libpmpi_count.so

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/cpmdc_pwscale.o: $(ROOT)/src/cpmdc_pwscale.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/cpmdc_pwscale_cuda_stub.o: $(ROOT)/src/cpmdc_pwscale_cuda_stub.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/pwscale_fft.o: $(ROOT)/bench/pwscale_fft.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/pwscale_harness.o: $(ROOT)/bench/pwscale_harness.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_pwscale_agreement.o: $(ROOT)/tests/test_pwscale_agreement.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/pmpi_count.o: $(ROOT)/bench/pmpi_count.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@

$(BUILD)/libpmpi_count.so: $(BUILD)/pmpi_count.o
	$(MPICC) -shared $^ -o $@

$(BUILD)/test_pwscale_collectives.o: $(ROOT)/tests/test_pwscale_collectives.c | $(BUILD)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_pwscale_collectives: $(BUILD)/test_pwscale_collectives.o $(BUILD)/pmpi_count.o $(BUILD)/cpmdc_pwscale.o $(BUILD)/cpmdc_pwscale_cuda_stub.o
	$(MPICC) $(LDFLAGS) $^ -lm -o $@

$(BUILD)/cpmdc-pwscale: $(BUILD)/pwscale_harness.o $(BUILD)/pwscale_fft.o $(BUILD)/cpmdc_pwscale.o $(BUILD)/cpmdc_pwscale_cuda_stub.o
	$(MPICC) $(LDFLAGS) $^ $(LIBS) -o $@

$(BUILD)/test_pwscale_agreement: $(BUILD)/test_pwscale_agreement.o $(BUILD)/cpmdc_pwscale.o $(BUILD)/cpmdc_pwscale_cuda_stub.o
	$(MPICC) $(LDFLAGS) $^ -lm -o $@

test: all
	$(BUILD)/test_pwscale_agreement
	mpirun --bind-to none --oversubscribe -np 2 $(BUILD)/test_pwscale_agreement
	mpirun --bind-to none --oversubscribe -np 4 $(BUILD)/test_pwscale_agreement
	mpirun --bind-to none --oversubscribe -np 8 $(BUILD)/test_pwscale_agreement
	mpirun --bind-to none --oversubscribe -np 2 $(BUILD)/test_pwscale_collectives
	mpirun --bind-to none --oversubscribe -np 4 $(BUILD)/test_pwscale_collectives

clean:
	rm -rf $(BUILD)
