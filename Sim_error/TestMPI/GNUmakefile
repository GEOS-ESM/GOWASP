#!/usr/bin/make 

BASELIBDIR = /discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.3_build1/x86_64-unknown-linux-gnu/ifort_13.1.2.183-mvapich2_1.8.1/Linux

_Ibufr	      = -I/discover/nobackup/rerrico/G5133/GEOSadas/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L/discover/nobackup/rerrico/G5133/GEOSadas/Linux/lib -lNCEP_bufr_r4i4

FC	= ifort

FCFLAGS	= -O3 -I. $(_Ibufr) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ibufr) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Ibufr)

mpLD	= mpif90
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lbufr)  

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  m_kinds.f90		        

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_RAD_PROFS	 = testMPI.F90
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


OBSFF90 = m_die.o 

ALL	= testMPI.x 
all: $(ALL)

# make excutable file for sim_rad.x
testMPI.x: $(LIB) $(OBSFF90) $(OBJS_RAD_PROFS) 
	$(mpLD) -o $@ $(OBJS_RAD_PROFS) $(OBSFF90) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x
lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

