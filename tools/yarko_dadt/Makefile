# makefile for swift_rmvsy and others yarko-integrators
# Miroslav Broz, miroslav.broz@email.cz
# Jul 3rd 2001

# compiler
f77 = gfortran

# options (profiling -pg, -mips4 -64 on SGI-HAL)
opt = -O3
FFLAGS = $(opt)
#FFLAGS = -64 -mips4 -O3
#LDFLAGS = -64 -mips4

# (-static on Linux, -non_shared for Digital Unix compiler)
static = -static
#static = -non_shared

# precompiled libraries
# libswift = libswift.a
libswift =

# lists of object files
obj = \
  io_init_spin.o \
  io_init_yarko.o \
  io_init_collision.o \
  yarko_omega.o \
  diurnal_dadt.o \
  seasonal_dadt.o \
  obliquity.o \
  io/io_init_param.o \
  io/io_init_pl.o \
  io/io_init_tp.o \
  io/io_getns.o \
  io/io_open.o \
  mvs/getacch/getacch_ir3.o \
  orbel/orbel_xv2el.o \
  util/util_exit.o \

inc = \
  const.inc \
  pr.inc \
  spin.inc \
  swift.inc \
  yarko.inc \
  io/io.inc

obj_all = $(obj)

# make rule
yarko_dadt : yarko_dadt.f $(obj) $(inc)
	$(f77) $(obj) $(opt) -o $@ $@.f $(libswift)

# comment out for SGI or Digital Unix
$(obj_all): %.o:%.f $(inc)
	$(f77) -c $(opt) -o $@ $<

clean : FORCE
	rm $(obj_all)
FORCE :

