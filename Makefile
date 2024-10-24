CC = gcc -g -ffreestanding#cc
CFLAGS = -O2 -fomit-frame-pointer -Wall -static
# CFLAGS = -Wall
# CFLAGS = -g -Wall -static
# CFLAGS = -pg -Wall -fprofile-arcs -static
LIBS = -lm 

PROG = PHYML
EXEC = phyml_multi 


DFLAG = -DUNIX -D$(PROG)

#rien : 
#	@echo Nothing to do ?...


######################################################################################################

OBJS = main.o utilities.o sigAnalyse.o filelocker.o optimiz.o ml.o bionj.o models.o free.o options.o simu.o eigen.o


$(EXEC) : $(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LIBS) $(CFLAGS)

clean :
	@rm *.o
######################################################################################################

eigen.o : eigen.c eigen.h
	$(CC) $(CFLAGS) $(DFLAG) -c eigen.c

simu.o : simu.c simu.h
	$(CC) $(CFLAGS) $(DFLAG) -c simu.c

ml.o : ml.c ml.h
	$(CC) $(CFLAGS) $(DFLAG) -c ml.c

utilities.o : utilities.c utilities.h
	$(CC) $(CFLAGS) $(DFLAG) -c utilities.c

sigAnalyse.o : sigAnalyse.c sigAnalyse.h
	$(CC) $(CFLAGS) $(DFLAG) -c sigAnalyse.c

filelocker.o : filelocker.c filelocker.h
	$(CC) $(CFLAGS) $(DFLAG) -c filelocker.c

optimiz.o : optimiz.c optimiz.h
	$(CC) $(CFLAGS) $(DFLAG) -c optimiz.c

bionj.o : bionj.c bionj.h
	$(CC) $(CFLAGS) $(DFLAG) -c bionj.c

main.o : main.c
	$(CC) $(CFLAGS) $(DFLAG) -c main.c

models.o : models.c models.h
	$(CC) $(CFLAGS) $(DFLAG) -c models.c

free.o : free.c free.h
	$(CC) $(CFLAGS) $(DFLAG) -c free.c

options.o : options.c options.h
	$(CC) $(CFLAGS) $(DFLAG) -c options.c

