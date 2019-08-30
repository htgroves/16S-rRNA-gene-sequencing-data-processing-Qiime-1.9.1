##### Codon Linux server

# Experiment number = HG0115

cd /data/lungen/16S_raw_data/Helen/HG0115 # go to folder containing sequencing data
mkdir A_sequence_and_mapping_files # make directory (folder) for the sequence and mapping files
cd A_sequence_and_mapping_files/ # go to new folder
gunzip *.gz # unzip files
module load qiime-1.9.0
extract_barcodes.py --input_type barcode_paired_end -f A_sequence_and_mapping_files/HG0115_S1_L001_I2_001.fastq -r A_sequence_and_mapping_files/HG0115_S1_L001_I2_001.fastq --bc1_len 8 --bc2_len 8 -o B_parsed_barcodes/ # Combines the forward and reverse barcodes that are in folder A and then sends them to folder B (B_parsed_barcodes)
trim_galore -a GGATTAGATACCCNNGTA -a2 CNCTTTANNCCCANT --length 150 --paired A_sequence_and_mapping_files/HG0115_S1_L001_R1_001.fastq A_sequence_and_mapping_files/HG0115_S1_L001_R2_001.fastq -o /data/lungen/16S_raw_data/Helen/HG0115/ # trimms the sequence reads
join_paired_ends.py -f HG0115_S1_L001_R1_001_val_1.fq -r HG0115_S1_L001_R2_001_val_2.fq -o /data/lungen/16S_raw_data/Helen/HG0115/C_joined_seqs/ -b /data/lungen/16S_raw_data/Helen/HG0115/B_parsed_barcodes/barcodes.fastq -j 200 -p 10  #joins together the forward and reverse reads where the minimum overlap length is 200bp and where <= 10% of the overlap can be misaligned 
mkdir run_quality_stats # make directy for quality stats to go in
nohup fastx_quality_stats -i C_joined_seqs/fastqjoin.join.fastq -o run_quality_stats/HG0115_quality_stats -Q33 # examins the read quality                                           
nohup split_libraries_fastq.py -i C_joined_seqs/fastqjoin.join.fastq -b C_joined_seqs/fastqjoin.join_barcodes.fastq -m A_sequence_and_mapping_files/HG0115_qiime_mapping_file_again.txt -o /data/lungen/16S_raw_data/Helen/HG0115/D_split_lib_all -q 29 --barcode_type 16 -r 10 -p 0.70 # demultiplexes the data in order to map each read back to its original sample. The read score must be above 30 ad if than more than 10 bases are lower than 30 the read is truncated. At least 70% of the combined read must be consecutively above 30 for the read to be retained
screen -R HG0115_phix # detaches screen
mkdir E_phix_removed # mkes directed for phix contaminantts and the no phix reads to go into
bwa aln -n 5 /data/lungen/microbiome_data/seq_process_scripts/phix_genome/phix_genome_tacked.fasta D_split_lib_all/seqs.fna > E_phix_removed/HG0115.sai # as the Phix amplicons have no barcode associated wit them, sometimes they take the barcodes of surrounding clusters. To get rid of this all the sequences are aligned with the Phix genome and any that match are removed.
bwa samse /data/lungen/microbiome_data/seq_process_scripts/phix_genome/phix_genome_tacked.fasta E_phix_removed/HG0115.sai D_split_lib_all/seqs.fna > phix_removed/mouse_gut_microbiome_bwa.sam # removing Phix
samtools view -F 4 -Sbh E_phix_removed/HG0115_bwa.sam > E_phix_removed/HG0115_phix_only.bam # removing Phix
bamtools convert -in E_phix_removed/HG0115_phix_only.bam -format fasta > E_phix_removed/HG0115_phix_only.fasta # removing Phix
samtools view -f 4 -Sbh E_phix_removed/HG0115_bwa.sam > E_phix_removed/HG0115_no_phix.bam # removing Phix
bamtools convert -in E_phix_removed/HG0115_no_phix.bam -format fasta > E_phix_removed/HG0115_no_phix.fasta # removing Phix
#ctrl AD to detach screen
nohup pick_open_reference_otus.py -i /data/lungen/16S_raw_data/Helen/HG0115/E_phix_removed/HG0115_no_phix.fasta -r /data/lungen/16S_reference_databases/silva_115/silva_115_database_final/sequences.fna -o /data/lungen/16S_raw_data/Helen/HG0115/F_uclust_open_ref_picked_otus_prefilter --prefilter_percent_id 0.6 -m uclust -a -O 8 -s 0.1 --suppress_taxonomy_assignment --suppress_align_and_tree # OTU picking using the open reference method uclust. Want to pick own reference OTUs and assign taxonomy using a different database so these arguments were suppressed.
pick_rep_set.py -i F_uclust_open_ref_picked_otus_prefilter/final_otu_map_mc2.txt -f E_phix_removed/HG0115_no_phix.fasta -o G_rep_set_uclust_prefilter/rep_set.fasta -m most_abundant -l G_rep_set_uclust_prefilter/log.txt # representative sequences were picked using the most abundant OTU. The mc2 OTU map argument removed any OTU with less than 2 sequences as these are likely to be sequencing errors
nohup align_seqs.py -i G_rep_set_uclust_prefilter/rep_set.fasta -t /data/lungen/16S_reference_databases/silva_115/silva_115_database_final/test_ref_align3.txt # sequences were aligned using Pynast
mv pynast_aligned H_pynast_aligned # aligned sequences moved into their own directory 
mkdir I_chimeric_seqs # new directiry created
nohup identify_chimeric_seqs.py -m ChimeraSlayer -i /data/lungen/16S_raw_data/Helen/HG0115/H_pynast_aligned/rep_set_aligned.fasta -a /data/lungen/16S_reference_databases/silva_115/silva_115_database_final/test_ref_align3.txt -o I_chimeric_seqs/chimeric_seqs.txt # chimera slayer identifies potential chimeric sequences. 
filter_fasta.py -f H_pynast_aligned/rep_set_aligned.fasta -o H_pynast_aligned/non_chimeric_rep_set_aligned.fasta -s I_chimeric_seqs/chimeric_seqs.txt -n # filters identified chimeric sequences out of aligned sequences
filter_alignment.py -i H_pynast_aligned/non_chimeric_rep_set_aligned.fasta -e 0.10 -g 0.80 -o H_pynast_aligned/ # the chimeria-free sequences were then screened for regions of high entropy as highly varibale regions are not very useful when peforming phylogenetic analysis
mkdir J_phylogenetic_tree # make directory for the tree to go in
nohup make_phylogeny.py -i H_pynast_aligned/non_chimeric_rep_set_aligned_pfiltered.fasta -o J_phylogenetic_tree/HG0115.tre # make a phylogenetic tree
nohup assign_taxonomy.py -i /data/lungen/16S_raw_data/Helen/HG0115/G_rep_set_uclust_prefilter/rep_set.fasta -t /data/lungen/16S_reference_databases/silva_115/silva_115_database_final/taxonomy.txt -r /data/lungen/16S_reference_databases/silva_115/silva_115_database_final/sequences.fna -o /data/lungen/16S_raw_data/Helen/HG0115/ -m rdp --rdp_max_memory 15000 # assign taxonomy using the RDP database. With large datasets the memory often has to be increased
mv rdp_assigned_taxonomy K_rdp_assigned_taxonomy # renames directory
mkdir L_otu_table
make_otu_table.py -i F_uclust_open_ref_picked_otus_prefilter/final_otu_map_mc2.txt -o L_otu_table/otu_table.biom -e I_chimeric_seqs/chimeric_seqs.txt -t K_rdp_assigned_taxonomy/rep_set_tax_assignments.txt # makes an OTU table with the taxonomy assignment but without the chimeric sequences
mkdir M_final_files_for_r # make the folder than will be exported into R for analysis
python /data/lungen/microbiome_data/seq_process_scripts/parse_rep_set_for_r.py G_rep_set_uclust_prefilter/rep_set.fasta M_final_files_for_r/rep_set_for_r.fasta # the headers of the representative sequences must be modified to be accepted into R and moved into the folder for R
cp J_phylogenetic_tree/HG0115.tre M_final_files_for_r # phylogenetic tree moved into the folder for R
biom convert -i L_otu_table/otu_table.biom -o otu_table_json.biom --table-type="OTU table" --to-json # OTU had to be given an attirubte to let R know what sort of OTU table it is.  Json is best accepted into R
cp otu_table_json.biom M_final_files_for_r # OTU table moved into folder for R

