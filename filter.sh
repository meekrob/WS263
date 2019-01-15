#!/usr/bin/env bash
# retain only "protein_coding" lines
# replace I-V,X with chrI,chrII,chrIII,chrIV,chrV,chrX
# replace MtDNA with chrM
# note: the sed commands have literal tabs in them (CTRL-V CTRL-I)
grep "protein_coding" /dev/stdin | sed 's/^I	/chrI	/; s/^II	/chrII	/; s/^III	/chrIII	/; s/^IV	/chrIV	/; s/^V	/chrV	/; s/^X	/chrX	/; s/^MtDNA	/chrM	/' 
