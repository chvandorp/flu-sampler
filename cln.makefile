NOTEBOOKS := $(wildcard scripts/*.ipynb)
CLEAN_NOTEBOOKS := $(NOTEBOOKS:.ipynb=.ipynb.cln)

all: $(CLEAN_NOTEBOOKS)

%.ipynb.cln: %.ipynb
	scripts/clear_output_ipynb.py $<
