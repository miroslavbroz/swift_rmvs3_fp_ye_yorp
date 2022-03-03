#!/usr/bin/awk -f

{
  gsub("cb_flt.inc",		"cb_flt_2ND.inc");
  gsub("cb_meanel.inc",		"cb_meanel_2ND.inc");
  gsub("cb_oscel.inc",		"cb_oscel_2ND.inc");
  gsub("cb_propel.inc",		"cb_propel_2ND.inc");
  gsub("proper.inc",		"proper_2ND.inc");

  gsub("/filter/",		"/filter_2ND/");
  gsub("/mean_elements/",	"/mean_elements_2ND/");
  gsub("/elements/",		"/elements_2ND/");
  gsub("/proper_elements/",	"/proper_elements_2ND/");
  gsub("/proper/",		"/proper_2ND/");

  gsub("filter_elmts_write",	"filter_elmts_write_2ND");
  gsub("fmft_call",		"fmft_call_2ND");
  gsub("fmft_gsout",		"fmft_gsout_2ND");
  gsub("io_init_filter",	"io_init_filter_2ND");
  gsub("io_init_proper",	"io_init_proper_2ND");
  gsub("io_write_filter",	"io_write_filter_2ND");
  gsub("io_write_fmft",		"io_write_fmft_2ND");
  gsub("io_write_fmft_propgs",	"io_write_fmft_propgs_2ND");
  gsub("io_write_proper",	"io_write_proper_2ND");
  gsub("proper_fmft",		"proper_fmft_2ND");
  gsub("proper_minmax2\\(",	"proper_minmax2_2ND(");
  gsub("proper_minmax\\(",	"proper_minmax_2ND(");
  gsub("proper_runavg",		"proper_runavg_2ND");
  gsub("proper_sigma\\(",	"proper_sigma_2ND(");
  gsub("proper_sigma2\\(",	"proper_sigma2_2ND(");
  gsub("proper_trojan",		"proper_trojan_2ND");
  gsub("proper_shift",		"proper_shift_2ND");

  print;
}

