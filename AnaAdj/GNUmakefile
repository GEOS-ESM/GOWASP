#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

FC	= ifort
FCFLAGS	= -O3 -I. $(_Inetcdf) 

mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lnetcdf) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  read_akbk.f90		\
		  m_prnt_test.f90	\
		  rw_nc4file.f90	


OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_CGRID	= enorm.f90
OBJS_CGRID	= $(addsuffix .o, $(basename $(SRCS_CGRID)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_CGRID) 

.SUFFIXES:
.SUFFIXES: .f90 .F90 .F .f .o

.f.o:
	$(FC) -c $(FCFLAGS) $<

.F.o:
	$(FC) -c $(FCFLAGS) $<

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.F90.o:
	$(mpLD) -c $(FCFLAGS) $<


ALL	= enorm.x 

all:	$(ALL)

# make excutable file for common_grid.x
enorm.x: $(LIB) $(OBJS_CGRID) 
	$(mpLD)   -o $@ $(OBJS_CGRID) $(_Lmine) $(LDFLAGS) 

clean:
	rm -f *.o *.mod *.x

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

