all: compile

compile:
	pdflatex pioneer-scan-v2.tex


dot2tex:
	./createDot.py

rundot: dot2tex compile

clean:
	rm pioneer-scan-v2.aux pioneer-scan-v2.bbl pioneer-scan-v2.blg

bibtex:
	bibtex pioneer-scan-v2

complete: dot2tex compile bibtex compile bibtex compile compile
	make compile
	make bibtex
	make compile

show:
	evince pioneer-scan-v2.pdf

edit:
	vim SMSdecomposition.tex

