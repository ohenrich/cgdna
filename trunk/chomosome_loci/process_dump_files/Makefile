
#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = c++

CFLAGS = -Wall -O3
INCLUDES =
LFLAGS = 
LIBS = 

# define the C source files
SRCS := $(wildcard *.cc)
OBJS = $(SRCS:.cc=.o)
MAIN = process_dumps.out

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo Done!!

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

capturec.o: capturec.h simFISH.h general.h
comparemaps.o: general.h comparemaps.h contactmap.h
contactmap.o: contactmap.h general.h
countbound.o: countbound.h general.h
general.o: general.h
pair_distance.o: pair_distance.h general.h
process_dumps.o: general.h process_dumps.h contactmap.h comparemaps.h rg.h
process_dumps.o: pair_distance.h countbound.h simFISH.h capturec.h
rg.o: rg.h general.h
simFISH.o: simFISH.h general.h
