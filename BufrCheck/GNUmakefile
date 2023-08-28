#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibufr	      = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L${GEOSPATH}/Linux/lib -lNCEP_bufr_r4i4

_Icrtm	      = -I${GEOSPATH}/src/NCEP_Shared/NCEP_crtm
_Lcrtm	      = -L${GEOSPATH}/Linux/lib -lNCEP_crtm

_Ibasic	      = -I../Lib_basic1
_Lbasic	      = -L../Lib_basic1 -lbasic1

_Iobsrw	      = -I../Lib_obsrw1
_Lobsrw	      = -L../Lib_obsrw1 -lobsrw1

FC	= ifort 
FCFLAGS	= -O3 -I.  $(_Iobsrw) $(_Ibasic) $(_Ibufr) $(_Icrtm)

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I.  $(_Iobsrw) $(_Ibasic) $(_Ibufr) $(_Icrtm)

#FCFLAGS = -g -O0 -check all -traceback -I.  $(_Iobsrw) $(_Ibasic) $(_Ibufr) $(_Icrtm)

mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		=  $(_Lobsrw) $(_Lbasic) $(_Lbufr) $(_Lcrtm) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  m_comp_stats.f90	\
                  m_count_types.f90	\
		  m_iasi_radiance.f90	\
		  m_read_loop.f90	\
                  m_saved_data.f90	
                  

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ERR	= check_bufr.f90

OBJS_ERR	= $(addsuffix .o, $(basename $(SRCS_ERR)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_ERR) 

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

ALL	= check_bufr.x 

all:	$(ALL)

# make excutable file for sim_rad.x
check_bufr.x: $(LIB) $(OBJS_ERR) 
	$(mpLD) -o $@ $(OBJS_ERR) $(_Lmine) $(LDFLAGS) 

clean:
	rm -f *.o *.mod *.x

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

