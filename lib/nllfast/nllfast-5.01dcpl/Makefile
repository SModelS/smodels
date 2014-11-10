all: compile chmod

tarball:
	tar czvf nllfast33.tar.gz nllfast-5.01dcpl.f *.grid Makefile

compile:
	gfortran -o nllfast_33TeV nllfast-5.01dcpl-33TeV.f

chmod:                                                                                
	touch gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out                        
	chmod 666 *.out

clean:
	rm -f nllfast_33TeV
