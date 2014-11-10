all: compile chmod

tarball:
	tar czvf nllfast13.tar.gz nllfast-3.0.f *.grid Makefile

compile:
	gfortran -o nllfast_13TeV nllfast-3.0-13TeV.f

chmod:                                                                                
	touch gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out                        
	chmod 666 *.out

clean:
	rm -f nllfast_13TeV
