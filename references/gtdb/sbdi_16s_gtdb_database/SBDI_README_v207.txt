# SBDI Sativa curated 16S GTDB database

## General information

Author: SBDI molecular data team
Contact e-mail: daniel.lundin@lnu.se, anders.andersson@scilifelab.se
DOI: 10.17044/scilifelab.14869077
License: CC BY 4.0
Version: R07-RS207-1
Categories: Bacteriology (310701), Microbial ecology (310703), Microbial genetics (310704), Medical bacteriology (320701), Medical microbiology not elsewhere classified (320799)
Item type: Dataset
Keywords: 16S rRNA, GTDB, DADA2, SBDI, Ampliseq
Funding: Curation of this data was funded by the Swedish Research Council (VR), grant number 2019-00242.

This readme file was last updated: 2021-09-02

Please cite as: Swedish Biodiversity Infrastructure (SBDI; 2021). SBDI Sativa curated 16S GTDB database. https://doi.org/10.17044/scilifelab.14869077

## Dataset description

The data in this [repository](https://doi.org/10.17044/scilifelab.14869077) is the result of vetting 16S sequences from the GTDB database release r207 (https://gtdb.ecogenomic.org/; Parks et al. 2018) with the Sativa program (Kozlov et al. 2016).

Files for the DADA2 (Callahan et al. 2016) methods `assignTaxonomy` and `addSpecies` are available, in three different versions each.
The `assignTaxonomy` files contain taxonomy for domain, phylum, class, order, family, genus and species.
(Note that it has been proposed that species assignment for short 16S sequences require 100% identity (Edgar 2018), so use species assignments with `assignTaxonomy` with caution.)
The versions differ in the maximum number of genomes that we included per species: 1, 5 or 20, indicated by "1genome", "5genomes" and "20genomes" in the file names respectively.
Using the version with 20 genomes per species should increase the chances to identify an exactly matching sequence by the `addSpecies` algorithm, while using a file with many genomes per species could potentially give biases in the taxonomic annotations at higher levels by `assignTaxonomy`.

There is also a fasta file with the original GTDB sequence names: gtdb-sbdi-sativa.r07rs207.20genomes.fna.gz.

All files are gzipped fasta files with 16S sequences, the assignTaxonomy associated with taxonomy hierarchies from domain to species whereas the `addSpecies` file have sequence identities and species names.

Taxonomical annotation of 16S amplicons using this data is available as an optional argument to the nf-core/ampliseq Nextflow workflow from version 2.1: `--dada_ref_taxonomy sbdi-gtdb` (https://nf-co.re/ampliseq; Straub et al. 2020).

The data will be updated circa yearly, after the GTDB database is updated.

### Curation

After download, sequences longer than 2000 basepairs and sequences containing undetermined bases ('N') were removed. 
Subsequently, sequences, as well as the reverse-complements of these, were aligned to the archaeal and bacterial SSU profiles from Barrnap (https://github.com/tseemann/barrnap) with hmmalign from HMMER (Eddy 2011). 
Sequences aligning to fewer than 1000 bases of their respective profile in both forward and reverse-complementary direction were deleted. 
For the sequences passing the above filters, the longest sequence in each genome was kept. 
For each species, a maximum of 20 sequences were selected, giving highest priority to sequences from GTDB species-representative genomes, and secondly longer sequences before shorter. 
These sequences were then analyzed with Sativa (Kozlov et al. 2016) and sequences misclassified at genus to phylum level were removed. 
For the remaining sequences, for each species, a maximum of 1, 5 and 20 sequences was selected, as before prioritizing sequences from GTDB species-representative genomes, and longer sequences before shorter. 
A Perl script for conducting filtering of sequences prior to and after Sativa analysis can be found in the `scripts` folder in the GitHub repo: https://github.com/biodiversitydata-se/sbdi-gtdb. 
Run `perl select_seq_sativa.pl --h` for documentation. 

## Version history

* v6 (20221007): Add missing fasta file with original GTDB names.

* v5 (20220902): Update README (this document)

* v4 (20220831): Update to GTDB R07-RS207 from R06-RS202

## References

Callahan, Benjamin J., Paul J. McMurdie, Michael J. Rosen, Andrew W. Han, Amy Jo A. Johnson, and Susan P. Holmes. 2016. “DADA2: High-Resolution Sample Inference from Illumina Amplicon Data.” Nature Methods 13 (7): 581–83. https://doi.org/10.1038/nmeth.3869.

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” PLoS Comput Biol 7 (10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195.

Edgar, Robert C. 2018. “Updating the 97% Identity Threshold for 16S Ribosomal RNA OTUs.” Bioinformatics 34 (14): 2371–75. https://doi.org/10.1093/bioinformatics/bty113.

Kozlov, Alexey M., Jiajie Zhang, Pelin Yilmaz, Frank Oliver Glöckner, and Alexandros Stamatakis. 2016. “Phylogeny-Aware Identification and Correction of Taxonomically Mislabeled Sequences.” Nucleic Acids Research 44 (11): 5022–33. https://doi.org/10.1093/nar/gkw396.

Parks, Donovan H., Maria Chuvochina, David W. Waite, Christian Rinke, Adam Skarshewski, Pierre-Alain Chaumeil, and Philip Hugenholtz. 2018. “A Standardized Bacterial Taxonomy Based on Genome Phylogeny Substantially Revises the Tree of Life.” Nature Biotechnology, August. https://doi.org/10.1038/nbt.4229.

Straub, Daniel, Nia Blackwell, Adrian Langarica-Fuentes, Alexander Peltzer, Sven Nahnsen, and Sara Kleindienst. 2020. “Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S RRNA (Gene) Amplicon Sequencing Pipeline.” Frontiers in Microbiology 11. https://doi.org/10.3389/fmicb.2020.550420.
