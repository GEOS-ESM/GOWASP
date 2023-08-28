#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

Lib1      = /home/rerrico/GOWASP_3/Lib_basic1
Lib2      = /home/rerrico/GOWASP_3/Lib_shmem1

_Icom1	      = -I$(Lib1)
_Lcom1	      = -L$(Lib1) -lbasic1

_Icom2	      = -I$(Lib2)
_Lcom2	      = -L$(Lib2) -lshmem1

FC	= ifort

FCFLAGS	= -O3 -I. $(_Inetcdf) $(_Icom1) $(_Icom2)

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Inetcdf) $(_Icom1) $(_Icom2) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Inetcdf) $(_Icom1) $(_Icom2)

mpLD	= mpif90
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lnetcdf) $(_Lcom1) $(_Lcom2)

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
                  m_write_ana_file.f90	

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_RAD_PROFS	 = NR2Egrid.F90
OBJS_RAD_PROFS	 = $(addsuffix .o, $(basename $(SRCS_RAD_PROFS)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_RAD_PROFS) 

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

ALL	= NR2Egrid.x 
all: $(ALL)

# make excutable file for sim_rad.x
NR2Egrid.x: $(LIB) $(OBJS_RAD_PROFS) 
	$(mpLD) -o $@ $(OBJS_RAD_PROFS) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a
lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

