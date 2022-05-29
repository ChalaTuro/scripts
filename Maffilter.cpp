## the input is sligned MAF file (*.maf) 
## the input was ordered to ref using maf_order.v10 in  MULTIZ tool. 
## maffilter param=option_file.cpp
####################################

DATA=PROJECTED2W11
input.file=$(DATA).maf
input.file.compression=none
input.format=Maf
output.log=$(DATA).filter.log          
maf.filter=     \
  Subset(        \
      species=(PtmSG1,Ptm16FRG073, Ptm17FRG089,PtmCad64,  PtmMur2,PttW11,Ptt17FRG026,Ptt9122,Ptt9139,Ptt9193,PttK103,PttNB29,PttNB85,Ptt773), \
      strict=yes, \
      keep=no,     \
      remove_duplicates=yes),      \
  MinBlockLength(min_length=100),   \ 
  MinBlockSize(min_size=3),          \
  RemoveEmptySequences(unresolved_as_gaps=yes), \   
  AlnFilter2(     \
    species=(PtmSG1,Ptm16FRG073, Ptm17FRG089,PtmCad64,  PtmMur2,PttW11,Ptt17FRG026,Ptt9122,Ptt9139,Ptt9193,PttK103,PttNB29,PttNB85,Ptt773), \
    window.size=10, \
    window.step=1,   	\
    max.gap=1,      		 \
    max.pos=1,       		\
    relative=no,      		 \
    missing_as_gap=yes, 	\
    file=check.block.AlnFilter2.aln), \
  Merge(       \
    species=(PtmSG1,Ptm16FRG073, Ptm17FRG089,PtmCad64, PtmMur2,PttW11,Ptt17FRG026,Ptt9122,Ptt9139,Ptt9193,PttK103,PttNB29,PttNB85,Ptt773), \
    dist_max=2, \   
    rename_chimeric_chromosomes=yes), \   
  Concatenate(minimum_size=10000000, ref_species=PttW11), \  
  SequenceStatistics(     \
    statistics=(                  \
      DiversityStatistics(    \
        ingroup=(PtmSG1,Ptm16FRG073, Ptm17FRG089,PtmCad64, PtmMur2,PttW11,Ptt17FRG026,Ptt9122,Ptt9139,Ptt9193,PttK103,PttNB29,PttNB85,Ptt773))), \
        ref_species=PttW11,     \ 
        file=$(DATA).diversity.stat.csv ), \
  OrderFilter(               \
    reference=PttW11,         \
    do_unsorted=discard,       \
    do_overlapping=discard),    \  
  Output(                                     \
    file=$(DATA).chromosome.ldhat.maf,         \
    compression=none,   \
    reference=PttW11,    \
    mask=yes),           \
  OutputAlignments(      \
    format=Fasta,        \
    mask=yes,             \
    coordinates=yes,      \
    compression=none,     \
    ldhat_header=yes,     \
    file=$(DATA).chromosome.ldhat.fasta), \          
  OutputAlignments(    \
    format=Fasta,     \
    mask=yes,         \
    coordinates=no,    \
    compression=none,   \
    ldhat_header=no,     \
    file=$(DATA).chromosome.no_ldhat.fasta), \  
  DistanceEstimation(  \
     method=ml,  \
     model=K80(kappa=2), \
     parameter_estimation=pairwise, \
     max_freq_gapps=0.3, \
     gaps_as_unresolved=yes, \
     profiler=profiles.txt, \
     message_handler=ml-message.txt, \
     extended_names=yes),  \
  DistanceBasedPhylogeny(  \
    method=bionj,  \
    dist_mat=MLDistance),  \
  OutputTrees(  \
    tree=BioNJ, \
    file=chromosomeLevelconcatenatedTree.dnd, \
    compression=none), \ 
  VcfOutput(                                    \
        file=chromLevelSNP.vcf,                   \
        reference=PttW11,                           \
        genotypes=(PtmSG1,Ptm16FRG073, Ptm17FRG089,PtmCad64,Ptt17FRG026,Ptt9122,Ptt9139,Ptt9193,PttK103,PttNB29,PttNB85,Ptt773), \        
        all=no)
#end of the analysis
