# Metisa_plana_target_mining

for blast can install from here: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

for uniprot database here: [https://www.uniprot.org/help/downloads](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)
 (make sure down all; both swiss and trembl also isoform)

             pipeline:
     
     quality check [00] -> clean file [01] -> blastp [02] -> threshold [03]
                           |
     contamination screening -> (bacteria etc) -> domain annotation ->  gene family/subfamily assignment
                           |
     candidate target selection -> Scoring candidates -> top-10 table -> biological justification
