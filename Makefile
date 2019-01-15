default: WS263.gtf

c_elegans.PRJNA13758.WS263.canonical_geneset.gtf.gz:
	wget ftp://ftp.wormbase.org/pub/wormbase/releases/current-production-release/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS263.canonical_geneset.gtf.gz

# this will be the filtered file
WS263.gtf: c_elegans.PRJNA13758.WS263.canonical_geneset.gtf.gz
	zcat $^ | ./filter.sh  > $@
