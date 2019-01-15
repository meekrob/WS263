#!/usr/bin/env python3
import gffutils,sys

# field to get
#ID_FIELD='transcript_id'
ID_FIELD='gene_id'

# open database
DB_PATH = '/projects/dcking@colostate.edu/support_data/annotations/wormbase/WS263.sqlite'
db = gffutils.FeatureDB(DB_PATH)

upstream = 2000
downstream = 2000

num_unmatched = 0
num_matched = 0

num_plus_strand_hits = 0
num_minus_strand_hits = 0
num_both_strand_hits = 0

strand_counts = {}

for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    line = line.strip()
    fields = line.split("\t")
    chrom, peak_s, peak_e = fields[0], fields[1], fields[2]
    wormbase_chrom = chrom.replace('chr','')

    ## Plus strand query ####################################################################################
    # Query a layout on the plus strand for a given transcript,peak pair                                    #
    # GIVEN: (start - upstream) < (start)                                                                   #
    # GIVEN:  peak_s < peak_e                                                                               #
    # EXPRESSION: peak_s < (start + downstream) AND peak_e > (start - upstream)                             #
    # ------------------------------------------------------------------------------------------------------#
    #           | Chromosomal orientation (+):                                                              #
    #    KEY    | Promoter: xxxxxxxx                                                                        #
    #           | Peak: ...........                                                                         #
    #           | Overlap:  *******                                                                         #
    # ------------------------------------------------------------------------------------------------------#
    #           :                        v-------------------promoter-----------------v                     #
    #           :                (start-upstream)xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx(start + down) >>>>>>>       #
    # Expression:                        ^--------------------------------------------^                     #
    #-----------: ------------------------------------------------------------------------------------------#
    #-----------: ----------------------------------------OVERLAPPING---------------------------------------#
    #   TRUE    :(peak_s ................******* peak_e)                                                    #
    #   TRUE    :                                           (peak_s *******************.... peak_e)         #
    #   TRUE    :(peak_s ................**********************************************..... peak_e)        #
    #   TRUE    :                           peak_s ********************* peak_e                             #
    #-----------: ------------------------------------NOT OVERLAPPING---------------------------------------#
    #  FALSE    :(peak_s ..... peak_e)                                                                      #
    #  FALSE    :                                                                      (peak_s ..... peak_e)#
    #########################################################################################################
    # EXPRESSION: peak_s - downstream < start  AND peak_e + upstream > start                                #
    #########################################################################################################
    feats = db.features_of_type("exon", (wormbase_chrom, int(peak_s) - downstream, int(peak_e) + upstream))
    #########################################################################################################
    plus_strand_hits = []
    for feat in feats:
        if feat.attributes['exon_number'][0] == '1':
            plus_strand_hits += feat.attributes[ ID_FIELD ]
            if feat.strand not in strand_counts: strand_counts[feat.strand] = 0
            strand_counts[feat.strand] += 1
            break

    ## Minus strand query ###################################################################################
    # query a layout on the plus strand for a given transcript,peak pair                                    #
    # Given: (start - upstream) < (start)                                                                   #
    # Given:  peak_s < peak_e                                                                               #
    # Expression: peak_s < (start+upstream) AND peak_e > (start-downstream)                                 #
    # --------------------------------------------------------------------------------                      #
    # Expression| Chromosomal orientation:                                                                  #
    # --------------------------------------------------------------------------------                      #
    # ------------------------------------------------------------------------------------------------------#
    #           | Chromosomal orientation (-):                                                              #
    #    KEY    | Promoter: xxxxxxxx                                                                        #
    #           | Peak: ...........                                                                         #
    #           | Overlap:  *******                                                                         #
    # ------------------------------------------------------------------------------------------------------#
    #           :                        v-------------------promoter-----------------v                     #
    #           :         <<<<<<< (start - down)xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx(start+upstream)          #
    # Expression:                        ^--------------------------------------------^                     #
    #-----------: ------------------------------------------------------------------------------------------#
    #-----------: ----------------------------------------OVERLAPPING---------------------------------------#
    #   TRUE    : peak_s ....................... peak_e                                                     #
    #   TRUE    :                                            peak_s ....................... peak_e          #
    #   TRUE    : peak_s ................................................................... peak_e         #
    #   TRUE    :                           peak_s ..................... peak_e                             #
    #-----------: ------------------------------------NOT OVERLAPPING---------------------------------------#
    #  FALSE    : peak_s ..... peak_e                                                                       #
    #  FALSE    :                                                                       peak_s ..... peak_e #
    #########################################################################################################
    # Expression: peak_s - upstream < start  AND peak_e + downstream > start                                #
    #########################################################################################################
    #feats = db.features_of_type("exon", (wormbase_chrom, int(peak_s) - upstream, int(peak_e) + downstream))
    #########################################################################################################
    #This needs to be the LAST exon

    minus_strand_hits = []

    if False:
        max_exon_number = -1
        exon_dict = {}
        for feat in feats:
            curr_exon_number = int(feat.attributes['exon_number'][0])
            max_exon_number = max(max_exon_number, curr_exon_number)
            exon_dict[curr_exon_number] = feat.attributes[ ID_FIELD ]

        if max_exon_number > -1:
            minus_strand_hits += exon_dict[max_exon_number]
            #print("- got hits", minus_strand_hits)

    if len(plus_strand_hits) > 0 and len(minus_strand_hits) > 0:
        num_both_strand_hits += 1
    elif len(plus_strand_hits) > 0:
        num_plus_strand_hits += 1
    else:
        num_minus_strand_hits += 1

    all_hits = plus_strand_hits + minus_strand_hits
    if len(all_hits) > 0:
        print(line, ",".join(all_hits), sep="\t")
        num_matched += 1
    else:
        print(line, "NA", sep="\t")
        num_unmatched += 1

    # Using file.flush() will cause buffered
    # output to write to the file but can 
    # negatively effect performance if done
    # frequently.
    if num_matched and num_matched % 250 == 0:
        sys.stdout.flush()
        print("strand counts", strand_counts, file=sys.stderr, flush=True)
    
print("Peaks matched: %d, unmatched: %d. %.2f fraction matched" % (num_matched, num_unmatched, num_matched / (num_matched+num_unmatched)), file=sys.stderr)
print("Plus strand hits: %d. Minus strand hits: %d. Hits on both strands: %d" % (num_plus_strand_hits, num_minus_strand_hits, num_both_strand_hits), file=sys.stderr)
print("strand counts", strand_counts, file=sys.stderr, flush=True)
