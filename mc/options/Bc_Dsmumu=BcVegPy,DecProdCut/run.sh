#!/bin/sh

SetupProject Gauss v45r5
gaudirun.py options.py $GAUSSOPTS/Gauss-2012.py $GAUSSOPTS/GenStandAlone.py $DECFILESROOT/options/14175001.py $LBBCVEGPYROOT/options/BcVegPy.py

