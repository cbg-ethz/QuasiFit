ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = quasifit

quasifit_SOURCES = quasifit.cpp quasifit.hpp
quasifit_CPPFLAGS = $(Eigen_CFLAGS) $(BOOST_CPPFLAGS)
quasifit_LDADD = $(GSL_LIBS) $(BOOST_CHRONO_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_SYSTEM_LIBS)

dist_doc_DATA = README.md
