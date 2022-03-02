#f77 = f77
f77 = gfortran
opt = -O3
obj = \
  io_read_hdr.o \
  io_read_hdr_r.o \
  io_read_hdr_r8.o \
  io_read_line.o \
  io_read_line_r.o \
  io_read_line_r8.o \

libswift = libswift.a
libswift =

follow2:  follow2.f $(obj)
	$(f77) $(opt) $(obj) -o $@ $@.f $(libswift)

$(obj): %.o:%.f
	$(f77) -c $(opt) $<

clean : FORCE
	rm $(obj)
FORCE :

