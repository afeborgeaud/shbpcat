# makefile for shbpcat
PROGS = mpi-shbpcat
FC	= mpif90
FC2=ifort
FFLAGS	= -O
#-heap-arrays
#debug
#FFLAGS = -O2 -heap-arrays -g -check all -fpe0 -warn -traceback -debug extended

SRC = calmat.f trialf.f glu2.f others.f dclisb.f dclisb3.f formpi.f mpi-shbpcat.f
SRC2 = calmat.f trialf.f others.f glu2.f dclisb.f dclisb3.f shbpcat.f
OBJS	= $(SRC:.f=.o)
OBJS2 =$(SRC2:.f=.o)
.SUFFIXES: .f .o

all:$(PROGS) shbpcat

mpi-shbpcat: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

shbpcat: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-shbpcat.o: mpi-shbpcat.f
	$(FC) $(FFLAGS) -c mpi-shbpcat.f -o $@


.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) mpi-shbpcat shbpcat workshbp
