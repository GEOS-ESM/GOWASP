#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GOWASP     = ${GOWASP_PATH}/Lib_basic1

_Ibasic	      = -I ${GOWASP}
_Lbasic	      = -L ${GOWASP} -lbasic1

FC	= ifort

FCFLAGS	= -O3 -I. $(_Ibasic) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ibasic) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Ibasic) 

mpLD	= mpif90
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lbasic) 
LDMATH	        = -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
                  sub_fit_cov.f90		\
                  sub_read_stats.f90

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_RAD_HCORR	 = eparams_rad_hcorr.f90
SRCS_CONV_HCORR	 = eparams_conv_hcorr.f90
SRCS_RAD_CHCORR	 = eparams_rad_chcorr.f90
SRCS_NEW_PARAMS	 = new_corr_err_params.f90
SRCS_GPS_VARS	 = gps_vars.f90

OBJS_RAD_HCORR	 = $(addsuffix .o, $(basename $(SRCS_RAD_HCORR)))
OBJS_CONV_HCORR	 = $(addsuffix .o, $(basename $(SRCS_CONV_HCORR)))
OBJS_RAD_CHCORR	 = $(addsuffix .o, $(basename $(SRCS_RAD_CHCORR)))
OBJS_NEW_PARAMS	 = $(addsuffix .o, $(basename $(SRCS_NEW_PARAMS)))
OBJS_GPS_VARS	 = $(addsuffix .o, $(basename $(SRCS_GPS_VARS)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_RAD_HCORR) $(SRCS_CONV_HCORR) $(SRCS_RAD_CHCORR) $(SRCS_NEW_PARAMS) $(SRCS_GPS_VARS)

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

ALL	= eparams_rad_hcorr.x eparams_conv_hcorr.x eparams_rad_chcorr.x new_corr_err_params.x gps_vars.x
all: $(ALL)

# make excutable file for sim_rad.x
eparams_rad_hcorr.x: $(LIB) $(OBJS_RAD_HCORR) 
	$(LD) -o $@ $(OBJS_RAD_HCORR) $(_Lmine) $(LDFLAGS)

eparams_conv_hcorr.x: $(LIB) $(OBJS_CONV_HCORR) 
	$(LD) -o $@ $(OBJS_CONV_HCORR) $(_Lmine) $(LDFLAGS)

eparams_rad_chcorr.x: $(LIB) $(OBJS_RAD_CHCORR) 
	$(LD) -o $@ $(OBJS_RAD_CHCORR) $(_Lmine) $(LDFLAGS)

new_corr_err_params.x: $(LIB) $(OBJS_NEW_PARAMS) 
	$(LD) -o $@ $(OBJS_NEW_PARAMS) $(_Lmine) $(LDFLAGS) $(LDMATH)

gps_vars.x:  $(LIB) $(OBJS_GPS_VARS) 
	$(LD) -o $@ $(OBJS_GPS_VARS) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a
lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

