.DELETE_ON_ERROR:
.SECONDARY:

BED := $(shell cat ../region.bed)
drug ?= rif

all: $(foreach region, $(BED), $(region).pickle.int)

%.pickle.int:
	./generate_pickle.py -b $* -i ~/iso 

