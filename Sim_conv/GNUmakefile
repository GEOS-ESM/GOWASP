#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibufr	      = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L${GEOSPATH}/Linux/lib -lNCEP_bufr_r4i4

_Ibasic	      = -I../Lib_basic1
_Lbasic	      = -L../Lib_basic1 -lbasic1

_Ishmem	      = -I../Lib_shmem1
_Lshmem	      = -L../Lib_shmem1 -lshmem1

_Iobsrw	      = -I../Lib_obsrw1
_Lobsrw	      = -L../Lib_obsrw1 -lobsrw1

FC	= ifort

FCFLAGS	= -O3 -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Inetcdf) $(_Ibufr) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Inetcdf) $(_Ibufr) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Inetcdf) $(_Ibufr) 

mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lobsrw) $(_Lshmem) $(_Lbasic) $(_Lnetcdf) $(_Lbufr) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
                  sonde_subs.f90		\
                  m_conv_types.f90		

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_CONV	 = create_conv.F90
OBJS_CONV	 = $(addsuffix .o, $(basename $(SRCS_CONV)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_CONV) 

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


OBSFF90 =  m_shmem_fields.o m_conv_obs.o

ALL	= create_conv.x 
all: $(ALL)

# make excutable file for sim_rad.x
create_conv.x: $(LIB) $(OBSFF90) $(OBJS_CONV) 
	$(mpLD) -o $@ $(OBJS_CONV) $(OBSFF90) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a
lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

