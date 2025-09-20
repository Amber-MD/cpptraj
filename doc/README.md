Cpptraj Documentation
=====================

The `GeneratePDFs.sh` script will attempt to build the Cpptraj documentation PDFs from the Lyx files. It requires md5sum and lyx. It checks the md5sum of each lyx file and if it differs from what is in `DocumentChecksums.txt` it will run Lyx to attempt to generate the PDF for that document only. 

If not using `GeneratePDFs.sh` you should ensure the \*.lyx file checksums are up to date in `DocumentChecksums.txt`.

Manual references are stored in `cpptraj.bib` in Bibtex format.
