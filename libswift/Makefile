# Makefile
# Makefile for libswift.a library.
# Miroslav Broz (miroslav.broz@email.cz), Mar 1st 2022

f77 = gfortran
cc = gcc

opt = -O3 -g -pg

obj = \
  anal/anal_energy.o \
  anal/anal_energy_write.o \
  anal/anal_jacobi.o \
  anal/anal_jacobi_write.o \
  bs/bs_der.o \
  bs/bs_int.o \
  bs/bs_step.o \
  coord/coord_b2h.o \
  coord/coord_b2h_tp.o \
  coord/coord_b2j.o \
  coord/coord_h2b.o \
  coord/coord_h2b_tp.o \
  coord/coord_h2j.o \
  coord/coord_j2b.o \
  coord/coord_j2h.o \
  coord/coord_vb2h.o \
  coord/coord_vb2h_tp.o \
  coord/coord_vh2b.o \
  coord/coord_vh2b_tp.o \
  coord/coord_vh2vj.o \
  discard/discard_massive.o \
  discard/discard_mass_merge.o \
  discard/discard_mass_peri.o \
  discard/discard_mass_reorder.o \
  discard/discard.o \
  discard/discard_peri.o \
  discard/discard_pl_close.o \
  discard/discard_pl.o \
  discard/discard_sun.o \
  io/io_discard_mass.o \
  io/io_discard_merge.o \
  io/io_discard_write.o \
  io/io_dump_param.o \
  io/io_dump_pl.o \
  io/io_dump_pl_symba.o \
  io/io_dump_tp.o \
  io/io_energy_write.o \
  io/io_getns.o \
  io/io_init_param.o \
  io/io_init_pl.o \
  io/io_init_pl_symba.o \
  io/io_init_tp.o \
  io/io_jacobi_write.o \
  io/io_lyap2_write.o \
  io/io_lyap_init.o \
  io/io_lyap_write.o \
  io/io_mvits_init.o \
  io/io_open.o \
  io/io_read_hdr.o \
  io/io_read_hdr_r.o \
  io/io_read_line.o \
  io/io_read_line_r.o \
  io/io_read_mass.o \
  io/io_splitname.o \
  io/io_write_frame.o \
  io/io_write_frame_r.o \
  io/io_write_hdr.o \
  io/io_write_hdr_r.o \
  io/io_write_line.o \
  io/io_write_line_r.o \
  io/io_write_mass.o \
  lyap2/lyap2_acc_tp.o \
  lyap2/lyap2_drift_dan.o \
  lyap2/lyap2_drift_kepu.o \
  lyap2/lyap2_drift_one.o \
  lyap2/lyap2_drift_tp.o \
  lyap2/lyap2_kepu_mikk.o \
  lyap2/lyap2_step.o \
  lyap2/lyap2_step_tp.o \
  lyap/lyap_renorm.o \
  lyap/lyap_step.o \
  lyap/lyap_step_sh.o \
  mvits/mvits_getacch_tp.o \
  mvits/mvits_init_rot.o \
  mvits/mvits_kickvh_tp.o \
  mvits/mvits_lstep_tp.o \
  mvits/mvits_rotate.o \
  mvits/mvits_rstep_tp.o \
  mvits/mvits_step.o \
  mvs/drift/drift_dan.o \
  mvs/drift/drift_kepmd.o \
  mvs/drift/drift_kepu_fchk.o \
  mvs/drift/drift_kepu_guess.o \
  mvs/drift/drift_kepu_lag.o \
  mvs/drift/drift_kepu_new.o \
  mvs/drift/drift_kepu.o \
  mvs/drift/drift_kepu_p3solve.o \
  mvs/drift/drift_kepu_stumpff.o \
  mvs/drift/drift.o \
  mvs/drift/drift_one.o \
  mvs/drift/drift_tp.o \
  mvs/getacch/getacch_ah0.o \
  mvs/getacch/getacch_ah1.o \
  mvs/getacch/getacch_ah2.o \
  mvs/getacch/getacch_ah3.o \
  mvs/getacch/getacch_ah3_tp.o \
  mvs/getacch/getacch_ir3.o \
  mvs/getacch/getacch.o \
  mvs/getacch/getacch_tp.o \
  mvs/kickvh/kickvh.o \
  mvs/kickvh/kickvh_tp.o \
  mvs/step/step_dkd.o \
  mvs/step/step_dkd_pl.o \
  mvs/step/step_dkd_tp.o \
  mvs/step/step_kdk.o \
  mvs/step/step_kdk_pl.o \
  mvs/step/step_kdk_tp.o \
  obl/obl_acc.o \
  obl/obl_acc_tp.o \
  obl/obl_pot.o \
  orbel/orbel_eget.o \
  orbel/orbel_ehie.o \
  orbel/orbel_ehybrid.o \
  orbel/orbel_el2xv.o \
  orbel/orbel_esolmd.o \
  orbel/orbel_fget.o \
  orbel/orbel_fhybrid.o \
  orbel/orbel_flon.o \
  orbel/orbel_scget.o \
  orbel/orbel_schget.o \
  orbel/orbel_xv2aeq.o \
  orbel/orbel_xv2el.o \
  orbel/orbel_zget.o \
  rmvs2/rmvs2_interp.o \
  rmvs2/rmvs2_interp_o.o \
  rmvs2/rmvs2_step.o \
  rmvs2/rmvs2_step_out.o \
  rmvs3/rmvs3_chk.o \
  rmvs3/rmvs3_elog.o \
  rmvs3/rmvs3_interp.o \
  rmvs3/rmvs3_step_in.o \
  rmvs3/rmvs3_step.o \
  rmvs3/rmvs3_step_out2.o \
  rmvs3/rmvs3_step_out.o \
  rmvs/rmvs_chk_ind.o \
  rmvs/rmvs_chk.o \
  rmvs/rmvs_elog.o \
  rmvs/rmvs_obl_acc.o \
  rmvs/rmvs_peri.o \
  rmvs/rmvs_step_in_mvpl.o \
  rmvs/rmvs_step_in.o \
  rmvs/rmvs_step_in_tp.o \
  rmvs/rmvs_step.o \
  rmvs/rmvs_step_out.o \
  tu4/tu4_getaccb.o \
  tu4/tu4_getaccb_tp.o \
  tu4/tu4_ldrift.o \
  tu4/tu4_ldrift_tp.o \
  tu4/tu4_step.o \
  tu4/tu4_vkickb.o \
  tu4/tu4_vkickb_tp.o \
  util/util_exit.o \
  util/util_hills1.o \
  util/util_hills.o \
  util/util_mass_peri.o \
  util/util_peri.o \
  util/util_version.o \

inc = \
  swift.inc \
  bs/bs.inc \
  io/io.inc \
  mvits/mvits.inc \
  rmvs/rmvs.inc \

lib = libswift.a

libswift : $(obj) $(inc)
	ar vr $(lib) $(obj)

$(obj) : %.o:%.f $(inc)
	$(f77) $(opt) -c -o $@ $<

clean : FORCE
	rm -f $(obj) $(lib)

FORCE :


