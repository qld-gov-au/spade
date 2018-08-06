# General options
OUTPUT_NAME := spade
OBJ_NAME := spade.obj
CC = cc
LINKER_FLAGS = -lgsl -lgslcblas -lm -lGL -lGLU -lglut

# Spade options
SPADE_CFLAGS = -g -std=c99 -Wall -Wextra -pedantic
SPADE_SOURCES = trackball.c

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
	$(CC) -c -g -std=c99 -Wall -Wextra -pedantic -o $(OBJ_NAME) spade.c
	$(CC) -g $(OBJ_NAME) $(OBJECTS) $(LINKER_FLAGS) -o $(OUTPUT_NAME)
	@echo "Build successful"

# Generate object files
$(SPADE_OBJECTS): %_SPADE.o: %.c
	$(CC) -c $(SPADE_CFLAGS) $< -o $@

$(MESCHACH_OBJECTS): %_MESCHACH.o: %.c
	$(CC) -c $(MESCHACH_CFLAGS) $< -o $@

# Regenerate executable
rebuild: clean build

# Remove build artifacts
clean:
	rm -f $(OUTPUT_NAME) $(SPADE_OBJECTS)


