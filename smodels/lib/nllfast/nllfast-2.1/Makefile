FCC=$(shell ../../check_fortran_compiler.sh)
PROGRAM=nllfast_8TeV

all: $(PROGRAM) gg.out sb.out sg.out ss.out st.out gdcpl.out sdcpl.out

$(PROGRAM):
	$(FCC) -o $(PROGRAM) nllfast-*.f
%.out:
	touch $@
	chmod 666 $@

clean:
	rm -f $(PROGRAM) *.out
