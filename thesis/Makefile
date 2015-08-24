all: build/thesis.pdf


TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build
                                                                                
build/thesis.pdf: tudothesis.cls *.tex references.bib | build
	test -f progress.log || echo 'time,words' > progress.log
	echo `date "+%Y-%m-%d %H:%M:%S"`,`detex thesis.tex | wc -w` >> progress.log
	python update_progress.py &
	texfot latexmk $(TeXOptions) thesis.tex


build:
	mkdir -p build/

clean:
	rm -rf build
