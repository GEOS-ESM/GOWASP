#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

FC	= ifort
CC	= cc
FCFLAGS	= -O3 -I. $(_Inetcdf) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Inetcdf)


mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LIB		= libshmem1.a
_Lmine		= -L. -lshmem1
_I		= -I.

SRCS_mylib	= \
		  MAPL_Sort.F90		        \
		  MAPL_ShmemMod.F90	        \
                  sort.c			\
		  m_die.F90		        

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ALL	= $(SRCS_mylib) 

.SUFFIXES:
.SUFFIXES: .f90 .F90 .F .f .c .o 

.c.o:
	$(CC) -c $<

.f.o:
	$(FC) -c $(FCFLAGS) $<

.F.o:
	$(FC) -c $(FCFLAGS) $<

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.F90.o:
	$(mpLD) -c $(FCFLAGS) $<

clean:	
	rm -f *.o *.mod *.x *.a

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

