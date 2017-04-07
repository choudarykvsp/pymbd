c_warnings_off = incompatible-pointer-types parentheses-equality shorten-64-to-32 \#warnings
f_warnings_off = maybe-uninitialized
CFLAGS = $(addprefix -Wno-,${c_warnings_off})
FFLAGS = $(addprefix -Wno-,${f_warnings_off})
FVENDOR = intelem #gnu95
CVENDOR = intelem #unix
FC = mpiifort #gfortran
## In case you are happy enough to work on the "gaia" cluster in Luxembourg you might need to also
## 'export LD_PRELOAD=</path/to/MKL/libraries/libmkl_core.so:/path/to/MKL/libraries/libmkl_intel_thread.so>'
## before running any calculation with this module
## (NOTE: also adapt libmkl_intel_thread.so in LD_PRELOAD according to threaded or sequential application)
LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl
