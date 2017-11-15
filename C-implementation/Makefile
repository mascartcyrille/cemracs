# COMPILER AND LINKER
CC				:=	gcc

# FLAGS
CFLAGS			:=	-Wall -Wextra -Werror -ansi -Wpedantic
LDFLAGS			:=
ADDFLAGS		:=	-lm
PROFFLAGS		:=	-pg
PARFLAGS		:=	-lpthread -fopenmp -Dlinux
OPTI			:=	-O3

# BINARY FOLDER AND TARGETS
BINFOLDER		:=	./bin/
SIM				:=	$(BINFOLDER)/stoch-pde
BTT				:=	$(BINFOLDER)/binToTxt
PROF			:=	$(BINFOLDER)/profiling
DEBUG			:=	$(BINFOLDER)/debug
PARALLEL		:=	$(BINFOLDER)/parallel
ALL				:=	$(SIM) $(BTT) $(PROF) $(DEBUG) $(PARALLEL)

# SOURCES
SRCFOLDER		:=	./src/
MAKEOUTFOLDER	:=	$(SRCFOLDER)/makeout/
BTTFOLDER		:=	$(SRCFOLDER)/binToTxt/
SIMFOLDER		:=	$(SRCFOLDER)/sim/
LIBFOLDER		:=	$(SRCFOLDER)/lib/
RNGDEP			:=	$(LIBFOLDER)/SFMT-src-1.5.1/SFMT.c \
					$(LIBFOLDER)/dSFMT-src-2.2.3/dSFMT.c
SIMDEP			:=	$(SIMFOLDER)/sim.c \
					$(RNGDEP)
SIMPARDEP		:=	$(SIMFOLDER)/sim_parallel_full.c \
					$(RNGDEP)

# $@: name of the target
# $<: name of first dependency
# $^: list of dependencies
# $?: list of dependencies more recent than target
# $*: name of the file (without suffix)

all: directories $(ALL)

$(PARALLEL): $(SIMPARDEP)
	$(CC) $(CFLAGS) $(OPTI) $(PARFLAGS) -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS)

$(PROF): $(SIMDEP)
	$(CC) $(CFLAGS) -fopenmp -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS) $(PROFFLAGS)	

$(DEBUG): $(SIMDEP)
	$(CC) $(CFLAGS) -fopenmp -p -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS)

$(SIM): $(SIMDEP)
	$(CC) $(CFLAGS) $(OPTI) -fopenmp -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS)

$(BTT): $(MAKEOUTFOLDER)/binToTxt.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(MAKEOUTFOLDER)/binToTxt.o: $(BTTFOLDER)/binToTxt.c
	$(CC) -o $@ -c $^ $(CFLAGS)
	
directories:
	@mkdir -p $(BINFOLDER)
	@mkdir -p $(MAKEOUTFOLDER)
	@mkdir -p $(BTTFOLDER)
	@mkdir -p $(SIMFOLDER)
	@mkdir -p $(LIBFOLDER)

clean:
	rm -rf *.o $(MAKEOUTFOLDER)/*.o

cleaner: clean
	rm -rf $(ALL)
