docs: CpptrajManual.pdf CpptrajDevelopmentGuide.pdf 

.lyx.pdf:
	lyx -batch --export pdf2 $<

CpptrajManual.pdf: CpptrajManual.lyx cpptraj.lyx
	lyx -batch --export pdf2 CpptrajManual.lyx

CpptrajDevelopmentGuide.pdf: CpptrajDevelopmentGuide.lyx
	lyx -batch --export pdf2 CpptrajDevelopmentGuide.lyx

clean:
	-/bin/rm -f CpptrajManual.pdf CpptrajDevelopmentGuide.pdf
