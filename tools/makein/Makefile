# Makefile for makein, baryc2 and laplac programs.
# Miroslav Broz (miroslav.broz@email.cz), Mar 13th 2003

#f77 = f77
f77 = gfortran

opt = -O3

obj_makein = pleph.o selcon.o a.o pv1.o
obj_makeinpl = a.o
obj_makeinsymba = a.o util/util_hills1.o
obj_laplac = rotat.o length.o
obj_lapinv = \
  rotat.o \
  orbel/orbel_el2xv.o \
  orbel/orbel_ehybrid.o \
  orbel/orbel_fhybrid.o \
  orbel/orbel_scget.o \
  orbel/orbel_schget.o \
  orbel/orbel_xv2el.o \
  orbel/orbel_zget.o \
  orbel/orbel_eget.o \
  orbel/orbel_ehie.o \
  orbel/orbel_esolmd.o \
  orbel/orbel_fget.o \
  orbel/orbel_flon.o \

obj_xvpl2el = a.o
obj_xvtp2el = a.o

obj = $(obj_makein) $(obj_makeinsymba) $(obj_laplac) $(obj_lapinv)
swift =

all : makein makeinpl makeinsymba baryc baryc2 laplac lapinv xvpl2el xvpl3el xvtp2el makeinceres

makein : makein.f $(obj_makein)
	$(f77) $(opt) $(obj_makein) -o $@ $<

makeinceres : makeinceres.f $(obj_makein)
	$(f77) $(opt) $(obj_makein) -o $@ $<

makeinpl : makeinpl.f $(obj_makeinpl)
	$(f77) $(opt) $(obj_makeinpl) -o $@ $<

makeinsymba : makeinsymba.f $(obj_makeinsymba)
	$(f77) $(opt) $(obj_makeinsymba) -o $@ $<

baryc : baryc.f
	$(f77) $(opt) -o $@ $<

baryc2 : baryc2.f
	$(f77) $(opt) -o $@ $<

laplac : laplac.f $(obj_laplac)
	$(f77) $(opt) $(obj_laplac) -o $@ $<

lapinv : lapinv.f $(obj_lapinv)
	$(f77) $(opt) $(obj_lapinv) -o $@ $< $(swift)

xvpl2el : xvpl2el.f $(obj_xvpl2el)
	$(f77) $(opt) $(obj_xvpl2el) -o $@ $<

xvpl3el : xvpl3el.f $(obj_xvpl2el)
	$(f77) $(opt) $(obj_xvpl2el) -o $@ $<

xvtp2el : xvtp2el.f $(obj_xvtp2el)
	$(f77) $(opt) $(obj_xvtp2el) -o $@ $<

$(obj) : %.o:%.f
	$(f77) $(opt) -c -o $@ $<

clean : FORCE
	rm $(obj)
FORCE :


