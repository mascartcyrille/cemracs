# COMPILER AND LINKER
CC				:=	gcc

# FLAGS
CFLAGS			:=	-Wall -Wextra -Werror -ansi -Wpedantic
LDFLAGS			:=
ADDFLAGS		:=	-lm
PROFFLAGS		:=	-pg
OPTI			:=	-O3

# BINARY FOLDER AND TARGETS
BINFOLDER		:=	./bin/
SIM				:=	$(BINFOLDER)/sim
BTT				:=	$(BINFOLDER)/binToTxt
PROF			:=	$(BINFOLDER)/profiling
DEBUG			:=	$(BINFOLDER)/debug
ALL				:=	$(SIM) $(BTT) $(PROF) $(DEBUG)

# SOURCES
SRCFOLDER		:=	./src/
MAKEOUTFOLDER	:=	$(SRCFOLDER)/makeout/
BTTFOLDER		:=	$(SRCFOLDER)/binToTxt/
SIMFOLDER		:=	$(SRCFOLDER)/sim/
LIBFOLDER		:=	$(SRCFOLDER)/lib/
SIMDEP			:=	$(SIMFOLDER)/sim.c \
					$(LIBFOLDER)/SFMT-src-1.5.1/SFMT.c \
					$(LIBFOLDER)/dSFMT-src-2.2.3/dSFMT.c

# $@: name of the target
# $<: name of first dependency
# $^: list of dependencies
# $?: list of dependencies more recent than target
# $*: name of the file (without suffix)

all: directories $(ALL)

$(PROF): $(SIMDEP)
	$(CC) $(CFLAGS) -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS) $(PROFFLAGS)	

$(DEBUG): $(SIMDEP)
	$(CC) $(CFLAGS) -p -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS)

$(SIM): $(SIMDEP)
	$(CC) $(CFLAGS) $(OPTI) -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ $(ADDFLAGS)

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
