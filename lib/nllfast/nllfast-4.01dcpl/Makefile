all: compile chmod

tarball:
	tar czvf nllfast14.tar.gz nllfast-4.01dcpl.f *.grid Makefile

compile:
	gfortran -o nllfast_14TeV nllfast-4.01dcpl-14TeV.f

chmod:                                                                                
	touch gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out                        
	chmod 666 *.out

clean:
	rm -rf nllfast_14TeV
