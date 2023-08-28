#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibasic	      = -I../Lib_basic1

FC	= mpiifort
FCFLAGS	= -O3 $(BIG_ENDIAN) -I. $(_Ibasic) $(_Inetcdf) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ibasic) $(_Inetcdf) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Ibasic) $(_Inetcdf) 

LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LIB		= libncdfods.a
_Lmine		= -L. -lncdfods
_I		= -I.

SRCS_mylib	= m_ods_RE.F90 	

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_ALL	= $(SRCS_mylib) 

.SUFFIXES:
.SUFFIXES: .F90 .o 

.F90.o:	
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

