all: compile chmod

tarball:
	tar czvf nllfast8.tar.gz nllfast-2.1.f *.grid Makefile

compile:
	gfortran -o nllfast_8TeV nllfast-2.1.f

chmod:
	touch gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out 
	chmod 666 *.out

clean:
	rm -f nllfast_8TeV
