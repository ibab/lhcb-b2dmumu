
green=$(shell echo -e "\033[0;32m")
reset=$(shell echo -e "\033[0;0m")

main.pdf: main.tex
	@latexmk -f -pdflatex='lualatex -halt-on-error' -pdf -g main 2>&1 | sed 's/Latexmk: /$(green)LATEXMK$(reset) /'
	@xdotool search --class mupdf key --window %@ r > /dev/null 2>&1


