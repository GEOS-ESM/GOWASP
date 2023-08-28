#!/usr/bin/make 
BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

_Ibufr	      = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L${GEOSPATH}/Linux/lib -lNCEP_bufr_r4i4

_Icrtm	      = -I${GEOSPATH}/src/NCEP_Shared/NCEP_crtm
_Lcrtm	      = -L${GEOSPATH}/Linux/lib -lNCEP_crtm

_Ibasic	      = -I../Lib_basic1
_Lbasic	      = -L../Lib_basic1 -lbasic1

_Ishmem	      = -I../Lib_shmem1
_Lshmem	      = -L../Lib_shmem1 -lshmem1

_Iobsrw	      = -I../Lib_obsrw1
_Lobsrw	      = -L../Lib_obsrw1 -lobsrw1

FC	= mpiifort 
FCFLAGS	= -O3 -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Ibufr) $(_Icrtm) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Ibufr) $(_Icrtm) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Ibufr) $(_Icrtm) 

mpLD	= mpiifort     #mpif90
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

mLib  = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
LDFLAGS		=  $(_Lobsrw) $(_Lshmem) $(_Lbasic) $(_Lbufr) $(_Lcrtm) $(mLib) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  m_obs_pert.f90	\
		  m_rf_diags_fields.f90	\
		  m_random_power.f90	\
		  m_rf_diags_power.f90	\
                  m_count_types.f90	\
		  m_rad2bt.f90		\
		  m_pert_loop.f90	\
                  m_shtrans_rf.f90	\
                  fft99.f90     	\
                  m_fft.f90     	\
                  gpsro_lat_func.f90  	\
                  m_random_gauss.f90	

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ERR	= create_error.F90

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

OBSFF90 = m_random_fields.o

ALL	= create_error.x 

all:	$(ALL)

# make excutable file for sim_rad.x
create_error.x: $(LIB) $(OBSFF90) $(OBJS_ERR) 
	$(mpLD) -o $@ $(OBJS_ERR) $(OBSFF90) $(_Lmine) $(LDFLAGS) 

clean:
	rm -f *.o *.mod *.x *.a

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

