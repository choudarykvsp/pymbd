FVENDOR = intelem #gnu95
CVENDOR = intelem #unix
FC = mpiifort
LDFLAGS = --link-lapack_opt
FFLAGS = -O3 #-check all
F2PY = f2py
BLDDIR = build
#AR = ar
#ARFLAGS = rcs
#IOLIBDIR = lib_io
-include system.mk

all: mbd

mbd: mbd.f90 \
	$(addprefix $(BLDDIR)/,mbd_interface.o mbd_helper.o)
	@mkdir -p $(BLDDIR)
	$(F2PY) -c --build-dir $(BLDDIR) --fcompiler=$(FVENDOR) \
	--f90exec=$(FC) --f90flags="$(FFLAGS)" --compiler=$(CVENDOR) \
	$(wordlist 2,3,$^) -m $@ $(LDFLAGS) $<
	rsync -a $(BLDDIR)/*.mod .

mbd_math: mbd_math.f90 $(addprefix $(BLDDIR)/,mbd.o mbd_interface.o mbd_helper.o)
	@mkdir -p $(BLDDIR)
	$(F2PY) -c --build-dir $(BLDDIR) --fcompiler=$(FVENDOR) \
		   --f90exec=$(FC) --f90flags="$(FFLAGS)" --compiler=$(CVENDOR) \
		   $(filter-out $<,$^) -m $@ $(LDFLAGS) $<
	rsync -a $(BLDDIR)/*.mod .

#mbd_scalapack: mbd_scalapack.f90 \
#	$(addprefix $(BLDDIR)/,mbd_interface.o mbd_helper.o) \
#	$(addprefix $(IOLIBDIR)/,scalapack_io.a)
#	@mkdir -p $(BLDDIR)
#	@mkdir -p $(IOLIBDIR)
#	$(F2PY) -c --build-dir $(BLDDIR) --fcompiler=$(FVENDOR) \
#	--f90exec=$(FC) --f90flags="$(FFLAGS)" --compiler=$(CVENDOR) \
#	$(wordlist 2,3,$^) $(IOLIBDIR)/scalapack_io.a -m $@ $(LDFLAGS) $<
#	rsync -a $(BLDDIR)/*.mod .

$(BLDDIR)/%.o: %.f90
	@mkdir -p $(BLDDIR)
	$(FC) -c -fPIC -J $(BLDDIR) $(FFLAGS) -o $@ $^

#$(IOLIBDIR)/scalapack_io.o: scalapack_io.f
#	@mkdir -p $(IOLIBDIR)
#	$(FC) -c -fPIC -o $@ $^

#$(IOLIBDIR)/scalapack_io.a: $(IOLIBDIR)/scalapack_io.o
#	@mkdir -p $(IOLIBDIR)
#	$(AR) $(ARFLAGS) $@ $^

clean:
	rm -f *.mod
	rm -rf $(BLDDIR)
#	rm -rf $(IOLIBDIR)

distclean: clean
	rm -f mbd.*so
	rm -rf mbd.*dSYM
