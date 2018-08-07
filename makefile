# Copyright 2016 State of Queensland
# This file is part of SPADE
# See spade.c, COPYING, COPYING.LESSER

# General options
OUTPUT_NAME := spade
CC = cc
LINKER_FLAGS = -lm -pthread

# Spade options
SPADE_CFLAGS = -g -lm -pthread -std=c99
SPADE_SOURCE_DIRS = initial machinery mathprop model optim util plotting
SPADE_SOURCES = \
	spade.c \
	common.c \
	parameters.c \
	arg.c \
	$(wildcard $(SPADE_SOURCE_DIRS:=/*.c)) \
	$(wildcard $(SPADE_SOURCE_DIRS:=/**/*.c)) \
	$(wildcard $(SPADE_SOURCE_DIRS:=/**/**/*.c))
SPADE_OBJECTS := $(patsubst %.c,%_SPADE.o,$(SPADE_SOURCES))

# Meschach options
MESCHACH_CFLAGS = -g
MESCHACH_SOURCES := \
	copy.c err.c matrixio.c memory.c vecop.c matop.c pxop.c submat.c \
	init.c otherio.c machine.c matlab.c ivecop.c version.c meminfo.c \
	memstat.c lufactor.c bkpfacto.c chfactor.c qrfactor.c solve.c \
	hsehldr.c givens.c update.c norm.c hessen.c symmeig.c schur.c \
	svd.c fft.c mfunc.c bdfactor.c
MESCHACH_SOURCES := $(addprefix meschach/, $(MESCHACH_SOURCES))
MESCHACH_OBJECTS := $(patsubst %.c,%_MESCHACH.o,$(MESCHACH_SOURCES)) 

# Computed options
OBJECTS := $(SPADE_OBJECTS) $(MESCHACH_OBJECTS)

# Default task
all: build

# Generate executable
build: clean $(OBJECTS)
	$(CC) $(OBJECTS) $(LINKER_FLAGS) -o $(OUTPUT_NAME)
	@echo "Build successful"

# Generate object files
$(SPADE_OBJECTS): %_SPADE.o: %.c
	$(CC) -c $(SPADE_CFLAGS) $< -o $@

$(MESCHACH_OBJECTS): %_MESCHACH.o: %.c
	$(CC) -c $(MESCHACH_CFLAGS) $< -o $@

# Regenerate executable
rebuild: clean build

# Run a sample project
test: rebuild
	./spade -fn karumba -alpha .09 -beta .09 -gamma 1.3 -iota-disabled .07 -kappa .1 -omega-disabled 160

# Remove build artifacts
clean:
	rm -f $(OUTPUT_NAME) $(OBJECTS)

