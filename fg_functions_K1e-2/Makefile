
obj = intep.o

intep_test: intep_test.f $(obj)
	gfortran $(obj) -o $@ $<

$(obj): %.o:%.f
	gfortran -c -o $@ $<

clean:
	rm $(obj)

