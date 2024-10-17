#######################################################################
#  Myprogram                                                          #
#######################################################################

OPT += -DTMAX=5

OPTIMIZE = -O2 -fopenmp

OPTIONS = $(OPTIMIZE) $(OPT)

EXEC = myprogram

OBJS = example.o

INCL = const.h example.h Makefile

CPPFLAGS = $(OPTIONS)

CC = g++

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS)

