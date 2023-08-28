#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

FC	= mpiifort
FCFLAGS	= -O3 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv 

#FCFLAGS = -g -O0 -check all -traceback 

mpLD	= mpiifort 
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LIB		= libbasic1.a
_Lmine		= -L. -lbasic1
_I		= -I.

SRCS_mylib	= \
		  compute_lon_stride.f90        \
		  find_name_subs.f90	        \
		  sub_zensun.f90	        \
		  interpolate_subs.f90	        \
                  set_file_name.f90             \
		  read_akbk.f90		        \
		  smooth_prof.f90	        \
		  m_kinds.f90		        \
		  m_parameters.f90	        \
		  m_nr_fields_info.f90	        \
		  m_sat_info_table.f90	        \
		  m_obs_error_table.f90	        \
                  m_time_compute.f90                  

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ALL	= $(SRCS_mylib) 

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

clean:	
	rm -f *.o *.mod *.x *.a



lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

