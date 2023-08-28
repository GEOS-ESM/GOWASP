#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

_Ibufr        = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4

_Ibasic	      = -I../Lib_basic1

FC	= mpiifort
FCFLAGS	= -O3 -I. $(_Ibasic) $(_Ibufr) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ibasic) $(_Ibufr) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Ibasic) $(_Ibufr) 

LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LIB		= libobsrw1.a
_Lmine		= -L. -lobsrw1
_I		= -I.

SRCS_mylib	= \
		  open_obs_file.f90 	\
		  m_conv_names.f90 	\
		  m_gpsro_names.f90	\
                  m_rad_obs_arrays.f90	\
                  m_bufr_rad.f90	\
                  m_bufr_conv.f90	\
		  m_bufr_gpsro.f90	

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ALL	= $(SRCS_mylib) 

.SUFFIXES:
.SUFFIXES: .f90 .f .o 

.f.o:
	$(FC) -c $(FCFLAGS) $<

.f90.o:	
	$(FC) -c $(FCFLAGS) $<

clean:	
	rm -f *.o *.mod *.x *.a

lib:	$(LIB) 

$(LIB): $(OBJS_mylib) 
	$(RM) $@
	$(AR) $@ $(OBJS_mylib) 

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

