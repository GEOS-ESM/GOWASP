#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibasic	      = -I../../Lib_basic1
_Lbasic	      = -L../../Lib_basic1 -lbasic1

_Incdfods     = -I../../Lib_ncdfods
_Lncdfods     = -L../../Lib_ncdfods -lncdfods

FCFLAGS	= -O3 -I. $(_Ibasic) $(_Incdfods) $(_Inetcdf) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ibasic) $(_Incdfods) $(_Inetcdf) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Ibasic) $(_Incdfods) $(_Inetcdf) 

FC	= ifort
mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lbasic) $(_Lncdfods) $(_Lnetcdf) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= 

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

APPSRCS := chcorr_get.f90       \
           chcorr_sumods.f90    \
           chcorr_tables.f90    \
           hcorr_sumods.f90     \
           hcorr_tables.f90     \
           prog1_corr.f90       \
           prog1_vcorr.f90      \
           vcorr_sumods.f90     \
           vcorr_tables.f90     \
           countobs_sum.f90     \
           countobs_tables.f90  \
           impacts_rad.f90      \
           impacts_conv.f90     \
           impacts_tables.f90   \
           read_j_value.f90     \
           jo_tables.f90

APPSRCS	:= $(wildcard $(APPSRCS))
APPOBJS	:= $(addsuffix .o, $(basename $(APPSRCS)))
APPS	:= $(addsuffix .x, $(basename $(APPSRCS)))

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

all : $(APPS)

%.x: $(LIB) %.o
	$(mpLD) -o $@ $*.o $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

