CC=gcc
CFLAGS=-Wall -Wextra -Werror -ansi -Wpedantic
LDFLAGS=
OPTI=-O3
EXEC=sim binToTxt

# $@: name of the target
# $<: name of first dependency
# $^: list of dependencies
# $?: list of dependencies more recent than target
# $*: name of the file (without suffix)

all: $(EXEC)

sim: sim.c SFMT-src-1.5.1/SFMT.c dSFMT-src-2.2.3/dSFMT.c
	$(CC) $(CFLAGS) $(OPTI) -DSFMT_MEXP=607 -DDSFMT_MEXP=521 -o $@ $^ -lm

binToTxt: binToTxt.o
	$(CC) -o $@ $^ $(LDFLAGS)

binToTxt.o: binToTxt.c
	$(CC) -o $@ -c $^ $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
