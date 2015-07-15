all: build/thesis.pdf


TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build
                                                                                
build/thesis.pdf: tudothesis.cls *.tex references.bib | build
	texfot latexmk $(TeXOptions) thesis.tex


build:
	mkdir -p build/

clean:
	rm -rf build
