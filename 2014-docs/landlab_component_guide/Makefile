.SUFFIXES:
.SUFFIXES: .tex .pdf

TEX_FILES = how_to_build_landlab_components.tex
PDF_FILES = ${TEX_FILES:.tex=.pdf}
TEMP_FILES = ${TEX_FILES:.tex=.log} \
					   ${TEX_FILES:.tex=.out} \
					   ${TEX_FILES:.tex=.aux}

# Use these commands
ifndef PDFLATEX
	PDFLATEX=pdflatex
endif
ECHO=/usr/bin/printf

# Use these compile flags
LATEXFLAGS=-halt-on-error -interaction=batchmode

PASS='\033[00;32mPASS\033[00m\n'
FAIL="\033[00;31mFAIL\033[00m\n"

.PHONY: clean

all: $(PDF_FILES)

.tex.pdf: $(TEX_FILES)
	@$(RM) -f $*.out 2> /dev/null || true
	@$(ECHO) "$(PDFLATEX) $?... "
	@($(PDFLATEX) $(LATEXFLAGS) $? > $*.out && $(ECHO) $(PASS)) || $(ECHO) $(FAIL) 
	@$(ECHO) "$(PDFLATEX) $?... "
	@($(PDFLATEX) $(LATEXFLAGS) $? >> $*.out && $(ECHO) $(PASS)) || $(ECHO) $(FAIL) 
	@$(ECHO) "$@ successfully created.\n"

clean:
	@$(RM) -f $(TEMP_FILES) $(PDF_FILES)

pdf: all clean
