#!/bin/csh -xf
setenv PATH_AB /archive/u/rerrico/G514osse/ana/Y2006/M07
cp $PATH_AB/*.ana.eta.20060*_00z.nc4 $NOBACKUP/ABfiles/.
cp $PATH_AB/*.ana.eta.20060*_12z.nc4 $NOBACKUP/ABfiles/.
cp $PATH_AB/*.bkg.eta.20060*_00z.nc4 $NOBACKUP/ABfiles/.
cp $PATH_AB/*.bkg.eta.20060*_12z.nc4 $NOBACKUP/ABfiles/.
