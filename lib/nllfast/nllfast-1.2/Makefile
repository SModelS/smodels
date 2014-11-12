PROGRAM=nllfast_7TeV

all: $(PROGRAM) gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out

$(PROGRAM):
	gfortran -o $(PROGRAM) nllfast-*.f
%.out:
	touch $@
	chmod 666 $@

clean:
	rm -f $(PROGRAM) *.out
