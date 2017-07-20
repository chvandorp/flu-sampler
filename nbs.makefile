CLEAN_NOTEBOOKS := $(wildcard scripts/*.ipynb.cln)
NOTEBOOKS := $(CLEAN_NOTEBOOKS:.ipynb.cln=.ipynb)

all: $(NOTEBOOKS)

%.ipynb: %.ipynb.cln
	cp $< $@
