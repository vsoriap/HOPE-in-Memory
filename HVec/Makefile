CFLAGS        = -g -march=native -std=c99
CFLAGS3_      = $(CFLAGS) -O3 #-fno-inline
CFLAGS0_      = $(CFLAGS) -O0
CFLAGS_      = $(CFLAGS) -O3
LDFLAGS      = -lm
LDFLAGS_     = $(LDFLAGS)
OMPFLAGS     = 

IINCL    = -I${EXTRAE_HOME}/include/
ILIBS    = -L$(EXTRAE_HOME)/lib 

PROGRAM_ = hirschberg
PROGS_   = $(PROGRAM_) $(PROGRAM_)_i

GCC_  = icc
GCC_  = gcc
VERSION_CODE ?= TEST

.PHONY: all 
all: $(PROGS_)



$(PROGRAM_): ./$(PROGRAM_).c
	$(GCC_) $(CFLAGS_) -fopenmp \
	$^ -o $@ $(LDFLAGS_)

$(PROGRAM_)_i: ./$(PROGRAM_).c
	$(GCC_) -D_EXTRAE_ $(CFLAGS_) $(OMPFLAGS) $(IINCL) \
	$^ -o $@ $(LDFLAGS_) $(ILIBS) -lseqtrace

clean:
	rm -f *.o $(PROGS_)
