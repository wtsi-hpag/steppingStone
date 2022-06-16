# Makefile for steppingStone 
CC= gcc
CFLAGS= -O4
LFLAGS= -lm 

SOURCES=step_fastq.c step_linkStones.c step_breakSort.c step_number.c step_breakProcess.c step_cleanProcess.c stepBreakPoint.c step_processStones.c step_shortReads.c step_edgeStones.c step_sortStones.c steppingStone.c

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,step-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) 


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe step-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp step-bin/steppingStone .
	@cp step-bin/stepPangenome .
	@cp step-bin/step_nation .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
