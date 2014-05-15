#!/bin/sh

if [[ "$1" == "--clean" ]]
then
	echo "Cleaning bootstrapped files"

	# ./
	rm -rf .deps/
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf autoscan-2.69.log
	rm -rf compile
	rm -rf config.guess
	rm -rf config.h
	rm -rf config.h.in
	rm -rf config.h.in~
	rm -rf config.log
	rm -rf config.status
	rm -rf config.sub
	rm -rf configure
	rm -rf configure.scan
	rm -rf depcomp
	rm -rf install-sh
	rm -rf libtool
	rm -rf ltmain.sh
	rm -rf m4/libtool.m4
	rm -rf m4/ltoptions.m4
	rm -rf m4/ltsugar.m4
	rm -rf m4/ltversion.m4
	rm -rf m4/lt~obsolete.m4
	rm -rf missing
	rm -rf stamp-h1

	# ./src
	rm -rf src/.deps/
	rm -rf src/Makefile
	rm -rf src/Makefile.in
	rm -rf src/quasifit
	rm -rf src/quasifit-quasifit.o
	rm -rf src/quasifit-io.o
	rm -rf src/quasifit-mcmc.o

	# Tarballs
	rm -rf quasifit-0.2.tar.bz2
	rm -rf quasifit-0.2.tar.gz

	# OS X files
	rm -rf .DS_Store src/.DS_Store
else
	echo "Bootstrapping Autotools"
	autoreconf -vif
fi
