#!/usr/bin/make

.PHONY: install_all Basic1_install Ncdfods_install Shmem1_install Obsrw1_install ComGrid_install Conv_install Error_install Gpsro_install Rad_install Satwind_install ErrProgs_install StatsProgs_install 
.PHONY: clean

install_all: Basic1_install Ncdfods_install Shmem1_install Obsrw1_install ComGrid_install Conv_install Error_install Gpsro_install Rad_install Satwind_install ErrProgs_install StatsProgs_install 

Basic1_install:
	$(MAKE) -C Lib_basic1 lib

Ncdfods_install:
	$(MAKE) -C Lib_ncdfods lib

Shmem1_install:
	$(MAKE) -C Lib_shmem1 lib

Obsrw1_install:
	$(MAKE) -C Lib_obsrw1 lib

ComGrid_install:
	$(MAKE) -C ComGrid all

Conv_install:
	$(MAKE) -C Sim_conv all

Error_install:
	$(MAKE) -C Sim_error all

Gpsro_install:
	$(MAKE) -C Sim_gpsro all

Rad_install:
	$(MAKE) -C Sim_rad all

Satwind_install:
	$(MAKE) -C Sim_satwind all

ErrProgs_install:
	$(MAKE) -C Tuning/ErrProgs all

StatsProgs_install:
	$(MAKE) -C Tuning/StatsProgs all

clean:
	$(MAKE) -C  Lib_basic1  $@
	$(MAKE) -C  Lib_ncdfods $@
	$(MAKE) -C  Lib_obsrw1  $@
	$(MAKE) -C  Lib_shmem1  $@
	$(MAKE) -C  BufrCheck   $@
	$(MAKE) -C  ComGrid     $@
	$(MAKE) -C  Sim_conv    $@
	$(MAKE) -C  Sim_error   $@
	$(MAKE) -C  Sim_gpsro   $@
	$(MAKE) -C  Sim_rad     $@
	$(MAKE) -C  Sim_satwind $@
	$(MAKE) -C  Tuning/ErrProgs $@
	$(MAKE) -C  Tuning/StatsProgs $@

