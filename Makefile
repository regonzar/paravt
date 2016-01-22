#       Parallel voronoi tessellation
#
#       CC       -- compiler (gcc)
#       CCOPTS1  -- compiler options (-O2)
#
#       commands: make & make clean

#######################################################################

#Default
CC     = mpicc
CCOPTS1 = -O2 -fpic

#######################################################################
# If MPI is located in a non default directory
INC   = 
LIBS  = 

CCOPTS2 = $(CCOPTS1)

# OBJS in execution frequency order.  CFILES after qhull.c are alphabetical
OBJS1 =  user_r.o global_r.o stat_r.o io_r.o geom2_r.o poly2_r.o \
       merge_r.o libqhull_r.o random_r.o geom_r.o poly_r.o qset_r.o mem_r.o userprintf_r.o usermem_r.o
OBJS = $(patsubst %.o,qhullsrc_r/%.o,$(OBJS1))$
OBJS2 = vars.o func.o iovt.o vt.o

all: paravt

.c.o:
	$(CC) -c $(CCOPTS1) $<

clean:
	rm -f *.o ../core paravt  libqhull.a \
            *.exe
	echo $(OBJS)

libqhull.a: $(OBJS)
	ar r libqhull.a $(OBJS1)
	-test -x /bin/ranlib -o -x /usr/bin/ranlib && ranlib libqhull.a

#qhull:  qhullsrc_r/unix_r.o libqhull.a
#	$(CC) -o qhull $(CCOPTS2) unix_r.o -L. -lqhull -lm 

paravt: paravt.c libqhull.a $(OBJS2)
	$(CC) paravt.c -o paravt $(CCOPTS2) $(INC) $(LIBS) $(OBJS2) -L. -I./qhullsrc_r/ -lqhull -lm 

# end of Makefile
                                                                                                            

