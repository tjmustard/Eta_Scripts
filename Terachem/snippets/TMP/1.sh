#!/bin/bash
for file in *.tcs
do
    #echo "YES!"
    base=`basename $file .tcs`
    head -n37 $file > $base.tmp
    echo "# Coordinate system for optimization algorithm (default: dlc)" >> $base.tmp
    echo "min_coordinates  dlc">>"$base".tmp
    echo "">>"$base".tmp
    echo "# Frequency of coordinates outputed">>"$base".tmp
    echo "min_dump 1">>"$base".tmp
    echo "">>"$base".tmp
    echo "# dump orbitals every MD step">>"$base".tmp
    echo "orbitalswrtfrq 1">>"$base".tmp
    echo "">>"$base".tmp
    echo "# Number of MD steps">>"$base".tmp
    echo "nstep 1000">>"$base".tmp
    echo "">>"$base".tmp
    echo "end">>"$base".tmp

    mv $base.tmp $file

done
