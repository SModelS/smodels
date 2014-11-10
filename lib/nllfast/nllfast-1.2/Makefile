all: compile chmod

tarball:
	tar czvf nllfast7.tar.gz nllfast-1.2.f *.grid Makefile

compile:
	gfortran -o nllfast_7TeV nllfast-1.2.f

chmod:                                                                                
	touch gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out                        
	chmod 666 *.out

clean:
	rm -f nllfast_7TeV
