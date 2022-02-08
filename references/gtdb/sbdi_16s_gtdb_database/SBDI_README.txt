# SBDI Sativa curated 16S GTDB database

## General information

Author: SBDI molecular data team
Contact e-mail: daniel.lundin@lnu.se
DOI: 10.17044/scilifelab.14869077
License: CC BY 4.0
Version: R06-RS202-1
Categories: Bacteriology (310701), Microbial ecology (310703), Microbial genetics (310704), Medical bacteriology (320701), Medical microbiology not elsewhere classified (320799)
Item type: Dataset
Keywords: 16S rRNA, GTDB, DADA2, SBDI, Ampliseq
Funding: Curation of this data was funded by the Swedish Research Council (VR), grant number 2019-00242.

This readme file was last updated: 2021-08-04

Please cite as: Swedish Biodiversity Infrastructure (SBDI; 2021). SBDI Sativa curated 16S GTDB database. https://doi.org/10.17044/scilifelab.14869077

## Dataset description

The data in this repository is the result of vetting 16S sequences from the GTDB database release R06-RS202 (https://gtdb.ecogenomic.org/; Parks et al. 2018) with the Sativa program (Kozlov et al. 2016).

Files for the DADA2 (Callahan et al. 2016) methods `assignTaxonomy` and `addSpecies` are available: gtdb-sbdi-sativa.r06rs202.assignTaxonomy.fna.gz and gtdb-sbdi-sativa.r06rs202.addSpecies.fna.gz.
There is also a fasta file with the original GTDB sequence names: gtdb-sbdi-sativa.r06rs202.fna.gz
All three files are gzipped fasta files with 16S sequences, the assignTaxonomy associated with taxonomy hierarchies from domain to genus whereas the addSpecies file have sequence identities and species names.

Taxonomical annotation of 16S amplicons using this data is available as an optional argument to the nf-core/ampliseq Nextflow workflow from version 2.1: `--dada_ref_taxonomy sbdi-gtdb` (https://nf-co.re/ampliseq; Straub et al. 2020).

The data will be updated circa yearly, when the GTDB database is updated.

### Curation

After download, sequences longer than 2000 basepairs and sequences containing undetermined bases ('N') were removed. 
Subsequently, sequences, as well as the reverse-complements of these, were aligned to the archaeal and bacterial SSU profiles from Barrnap (https://github.com/tseemann/barrnap) with hmmalign from HMMER (Eddy 2011). 
Sequences aligning to fewer than 1000 bases of their respective profile in both forward and reverse-complementary direction were deleted. 
For the sequences passing the above filters, the longest sequence in each genome was kept. 
For each species, a maximum of 5 sequences was selected, prioritizing sequences from GTDB species-representative genomes, and longer sequences before shorter. 
The remaining sequences were then analyzed with Sativa (Kozlov et al. 2016) and sequences misclassified at genus to phylum level were removed.
A Perl script for conducting filtering of sequences prior to and after Sativa analysis can be found in the `scripts` folder in the GitHub repo: https://github.com/biodiversitydata-se/sbdi-gtdb. 
Run `perl select_seq_sativa.pl --h` for documentation. 

## References

Callahan, Benjamin J., Paul J. McMurdie, Michael J. Rosen, Andrew W. Han, Amy Jo A. Johnson, and Susan P. Holmes. 2016. “DADA2: High-Resolution Sample Inference from Illumina Amplicon Data.” Nature Methods 13 (7): 581–83. https://doi.org/10.1038/nmeth.3869.

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” PLoS Comput Biol 7 (10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195.

Kozlov, Alexey M., Jiajie Zhang, Pelin Yilmaz, Frank Oliver Glöckner, and Alexandros Stamatakis. 2016. “Phylogeny-Aware Identification and Correction of Taxonomically Mislabeled Sequences.” Nucleic Acids Research 44 (11): 5022–33. https://doi.org/10.1093/nar/gkw396.

Parks, Donovan H., Maria Chuvochina, David W. Waite, Christian Rinke, Adam Skarshewski, Pierre-Alain Chaumeil, and Philip Hugenholtz. 2018. “A Standardized Bacterial Taxonomy Based on Genome Phylogeny Substantially Revises the Tree of Life.” Nature Biotechnology, August. https://doi.org/10.1038/nbt.4229.

Straub, Daniel, Nia Blackwell, Adrian Langarica-Fuentes, Alexander Peltzer, Sven Nahnsen, and Sara Kleindienst. 2020. “Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S RRNA (Gene) Amplicon Sequencing Pipeline.” Frontiers in Microbiology 11. https://doi.org/10.3389/fmicb.2020.550420.
