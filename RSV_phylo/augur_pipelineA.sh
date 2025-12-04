## augur pipeline
## either use nextstrain shell envi
set -e
subtype="A"
localFasta="data/All_RSV_${subtype}_genomes_aln_qc_masked1.fasta"
#contextFasta="data/gisaid_wBlast_2000_2023_countryQuota100_${subtype}_E20.fasta"
contextFasta="../GISAID_RSV_A/gisaid_rsv_2025_full_EPI1.fasta"

ref="KT992094.1"
refgb="data/${ref}.gb"
version=3

seqs="data/unalignment${subtype}_${version}.fasta"
#seqs="data/gisaid_full_local_${subtype}_unaligned.fasta"
ali="data/alignment${subtype}_${version}.fasta"
tree="tree/tree_${subtype}_${version}.nwk"

meta="data/meta_${subtype}c.csv"
timetree="results/timetree_${subtype}_${version}.nwk"

branches=results/branch_lengths_${subtype}_${version}.json
traits=results/traits_${subtype}_${version}.json
muts=results/nt_muts_${subtype}_${version}.json
aamuts=results/aa_muts_${subtype}_${version}.json

output=auspice/rsv_${subtype}_${version}c.json
#cols="Lineage Nationality ARI.SARI Gender Results Referring.Facility"
#cat $localFasta $contextFasta > $seqs
## alternatively:
#cat data/unalignedA_local.fasta ../GISAID_RSV_A/gisaid_rsv_2025_full.fasta > data/gisaid_full_local_A_unaligned.fasta

#echo "################### ALIGNMENT #######################"
#augur align --sequences $seqs  --nthreads 30 --method mafft --reference-name $ref --output $ali
#translates to
#mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 32 alignment.fasta.to_align.fasta 1> alignment.fasta 2> alignment.fasta.log


## from root directory
#echo "################### TREE Build #######################"
#augur tree --alignment $ali --nthreads 16 --output $tree

echo "################### TREE Refinement #######################"
augur refine --alignment $ali --tree $tree --metadata $meta --output-tree $timetree --output-node-data $branches --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 5
#--clock-filter-iqd 100 kicks out ref KT992094.1 ?!?!?

echo "################### TIME-TREE Build (traits) #######################"
augur traits --tree $timetree --metadata $meta --output-node-data $traits --columns country --confidence

echo "################### ANCESTRAL #######################"
augur ancestral --tree $timetree --alignment $ali --output-node-data $muts --inference joint

echo "################### TRANSLATION #######################"
augur translate --tree $timetree --ancestral-sequences $muts --reference-sequence $refgb --output-node-data $aamuts

echo "################### FREQUENCIES #######################"
augur frequencies --method kde --metadata $meta --min-date 2023-01-01 --max-date 2024-12-31 --tree $timetree \
      --include-internal-nodes True --output results/cladefreq.json
      
echo "################### EXPORT to JSON  #######################"
augur export v2 --title "RSV ${subtype}" --maintainers "Henschel, Senghore, Zhou, Everett" --tree $timetree --metadata $meta --node-data $branches $traits $muts $aamuts \
      --colors config/colors.tsv --auspice-config config/auspice_config.json \
      --lat-longs config/lat_long_all.tsv \
      --color-by country Gender Nationality Facility ARI.SARI AgeGroup clade --output $output
