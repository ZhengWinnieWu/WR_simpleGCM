include intel.mk

wr_hs.x: libfms_spectral_solo.a
        $(LD) $^ $(LDFLAGS) -o $@

libfms_spectral_solo.a:  FORCE
        make  NETCDF=3  -f Makefile.fms_spectral_solo $@

FORCE:

clean:
        make  NETCDF=3  -f Makefile.fms_spectral_solo clean

localize:
        make  NETCDF=3  -f Makefile.fms_spectral_solo localize
