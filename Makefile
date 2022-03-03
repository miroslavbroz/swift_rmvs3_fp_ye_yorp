# Makefile
# Makefile for swift_rmvs3_fp_ye_yorp integrator.
# Miroslav Broz (miroslav.broz@email.cz), Feb 14th 2008

f77 = gfortran
cc = gcc

opt = -O3

lib = -lswift -L.

obj = \
  io/io_close.o \
  io/io_write_frame_r8_IO.o \
  io/io_write_hdr_r8.o \
  io/io_write_line_r8.o \
  filter/bessi0.o \
  filter/dft.o \
  filter/lsm.o \
  filter/filter_create.o \
  filter/filter_elmts.o \
  filter/filter_elmts_write.o \
  filter/filter_filter.o \
  filter/filter_shift.o \
  filter/fmft_call.o \
  filter/fmft_gsout.o \
  filter/io_init_filter.o \
  filter/io_init_proper.o \
  filter/io_write_filter.o \
  filter/io_write_proper.o \
  filter/io_write_fmft.o \
  filter/io_write_fmft_propgs.o \
  filter/proper_fmft.o \
  filter/proper_minmax.o \
  filter/proper_minmax2.o \
  filter/proper_runavg.o \
  filter/proper_shift.o \
  filter/proper_sigma.o \
  filter/proper_sigma2.o \
  filter/proper_trojan.o \
  filter2/filter_elmts_write_2ND.o \
  filter2/fmft_call_2ND.o \
  filter2/fmft_gsout_2ND.o \
  filter2/io_init_filter_2ND.o \
  filter2/io_init_proper_2ND.o \
  filter2/io_write_filter_2ND.o \
  filter2/io_write_fmft_2ND.o \
  filter2/io_write_fmft_propgs_2ND.o \
  filter2/io_write_proper_2ND.o \
  filter2/proper_fmft_2ND.o \
  filter2/proper_minmax2_2ND.o \
  filter2/proper_minmax_2ND.o \
  filter2/proper_runavg_2ND.o \
  filter2/proper_shift_2ND.o \
  filter2/proper_sigma2_2ND.o \
  filter2/proper_sigma_2ND.o \
  filter2/proper_trojan_2ND.o \
  misc/ang_eq.o \
  misc/ang_sgndot.o \
  misc/arr_avg.o \
  misc/arr_maxind.o \
  misc/arr_minind.o \
  misc/arr_minmax.o \
  misc/arr_shift.o \
  misc/dot_product3.o \
  misc/intep.o \
  misc/isgn.o \
  misc/length.o \
  misc/normalize.o \
  misc/ran1.o \
  misc/read_f.o \
  misc/read_i.o \
  misc/read_l.o \
  misc/read_s.o \
  misc/rotate_around_axis.o \
  misc/vector_product.o \
  misc/zero2pi.o \
  mvs/drift/drift.o \
  mvs/drift/drift_tp.o \
  mvs/getacch/getacch_ah3_tp.o \
  mvs/getacch/getacch_tp.o \
  mvs2/step_kdk2.o \
  mvs2/step_kdk2_pl.o \
  mvs2/step_kdk2_tp.o \
  rmvs3/rmvs3_chk.o \
  rmvs3/rmvs3_interp.o \
  rmvs3/rmvs3_step.o \
  rmvs3/rmvs3_step_out.o \
  util/util_version.o \
  yarko/discard_meana.o \
  yarko/disrupt.o \
  yarko/disrupt_write.o \
  yarko/getacc_yarko.o \
  yarko/io_dump_spin.o \
  yarko/io_init_collision.o \
  yarko/io_init_spin.o \
  yarko/io_init_yarko.o \
  yarko/io_init_yorp.o \
  yarko/omega_crit.o \
  yarko/reorient.o \
  yarko/reorient_write.o \
  yarko/yarko_omega.o \
  yarko/yarko_seasonal.o \
  yarko/yorp_evolve.o \

objc = \
  filter/fmft.o \
  filter/nrutil.o \

inc = \
  swift.inc \
  version.inc \
  filter/filter.inc \
  filter/proper.inc \
  filter/cb_flt.inc \
  filter/cb_meanel.inc \
  filter/cb_oscel.inc \
  filter/cb_propel.inc \
  filter2/proper_2ND.inc \
  filter2/cb_flt_2ND.inc \
  filter2/cb_meanel_2ND.inc \
  filter2/cb_oscel_2ND.inc \
  filter2/cb_propel_2ND.inc \
  yarko/spin.inc \
  yarko/yarko.inc \
  yarko/yorp.inc \

all: main/swift_rmvs3_fp_ye_yorp main/swift_mvs2_fp_ye_yorp main/swift_mvs2_fp2

main/swift_rmvs3_fp_ye_yorp: main/swift_rmvs3_fp_ye_yorp.f libswift.a $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(lib)

main/swift_mvs2_fp_ye_yorp: main/swift_mvs2_fp_ye_yorp.f libswift.a $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(lib)

main/swift_mvs2_fp2: main/swift_mvs2_fp2.f libswift.a $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(lib)

$(obj) : %.o:%.f $(inc)
	$(f77) $(opt) -c -o $@ $<

$(objc) : %.o:%.c
	$(cc) $(opt) -c -o $@ $<

libswift.a :
	cd libswift; make

clean : FORCE
	rm -f $(obj) $(objc)

FORCE :


