ml Python

for ((n=0;n<10000;n++))
do
  python3 Randomize_GO_terms_isoseq_Noremoval.py Ensembl.Genes.GO.Terms.txt 6844 Random.GO.Terms.Novel.Isoforms.with.GO.terms.Distribution.txt ${n}
  echo ${n}
done
