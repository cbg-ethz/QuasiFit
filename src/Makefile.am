PRECISION_LIBS = @PRECISION_LIBS@

bin_PROGRAMS = quasifit

quasifit_SOURCES = io.cpp mcmc.cpp quasifit.cpp quasifit.hpp ezETAProgressBar.hpp console_buf.hpp vars.hpp
quasifit_CXXFLAGS = $(EIGEN_CFLAGS) $(BOOST_CPPFLAGS)
quasifit_LDADD = $(GSL_LIBS) $(BOOST_CHRONO_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_SYSTEM_LIBS) $(PRECISION_LIBS)
