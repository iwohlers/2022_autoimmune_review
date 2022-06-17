
EFO_IDS = []
with open("data/autoimmune_efo_ids.tsv","r") as f_in:
    for line in f_in:
        EFO_IDS.append(line.strip())


rule all:
    input: "results/autoimmune_pgs_all_metadata_scores.csv",
           "results/autoimmune_pgs_all_metadata_evaluation_sample_sets.csv",
           "results/autoimmune_pgs_all_metadata_performance_metrics.csv",
           "results/autoimmune_pgs_all_metadata_score_development_samples.csv",
           "results/autoimmune_pgs_all_metadata_evaluation_sample_sets.csv",
           "results/autoimmune_pgs_count.tsv",
           expand("results/{efo_id}_studycounts_study.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_study.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_snps.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_reportedgenes.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_snpidcurrent.tsv",efo_id=EFO_IDS),
           "results/autoimmune_UK_Biobank_master_file.tsv",

# To run this, "results/autoimmune_pgs_all_metadata_scores.csv" needs to be 
# edited so text within "" doesn't contain a ","; the file should be named
# results/autoimmune_pgs_all_metadata_scores.csv
rule pgs_disease_specific_all:
    input:
           expand("results/{efo_id}_pgsids.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgs_all_metadata_score_development_samples.csv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgs_all_metadata_performance_metrics.csv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgs_all_metadata_evaluation_sample_sets.csv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgs_performance_metrics_count.csv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgseval_pubcount.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgs_pubcount.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_pgscombined_pubcount.tsv",efo_id=EFO_IDS)
           
 
rule get_header:
    input: "data/gwas_catalog/2022/05/23/gwas-catalog-associations_ontology-annotated.tsv"
    output: "results/header_gwas-catalog-associations_ontology-annotated.tsv"
    shell: "head -n 1 {input} > {output}"

# We select all associations that have only one mapped trait URI (EFO ID)
# and for which this one EFO ID is within the list of manually obtained 
# Autoimmune disease EFO IDs
# Results in 5,022 associations
rule select_autoimmune_disease_by_efo:
    input: "data/autoimmune_efo_ids.tsv",
           "data/gwas_catalog/2022/05/23/gwas-catalog-associations_ontology-annotated.tsv"
    output: "results/autoimmune_gwas-catalog-associations_ontology-annotated.tsv"
    run:
        efo_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                efo_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split("\t")
                 mapped_trait = s[35]
                 if (not ',' in mapped_trait and mapped_trait.split("/")[-1] in efo_ids) or line[:10] == "DATE ADDED":
                    f_out.write(line)

# Count the number of associations per trait (based on EFO ID)
rule count_associations_per_disease:
    input: "results/autoimmune_gwas-catalog-associations_ontology-annotated.tsv"
    output: "results/autoimmune_association_count.tsv"
    shell: "cat {input} | cut -f 36 | sort | uniq -c | sort -k1 -n -r > {output}"

# Get disease-specific lists
rule get_disease_specific_associations:
    input: "results/autoimmune_gwas-catalog-associations_ontology-annotated.tsv"
    output: "results/{efo_id}-associations_ontology-annotated.tsv"
    run: 
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split("\t")
                 mapped_trait = s[35]
                 if mapped_trait.split("/")[-1] == wildcards.efo_id or line[:10] == "DATE ADDED":
                    f_out.write(line)    
    
rule get_all_efo_ids:
    input: expand("results/{efo_id}-associations_ontology-annotated.tsv",efo_id=EFO_IDS)

# Count the number of GWAS catalog study accessions (column 37, 
# "STUDY ACCESSION") per trait (based on EFO ID)
rule count_studies_disease_specific:
    input: "results/{efo_id}-associations_ontology-annotated.tsv"
    output: "results/{efo_id}_counts_study.tsv"
    shell: "cat {input} | cut -f 37 | sort | uniq -c | sort -k1 -n -r > {output}"

# Count the number of GWAS catalog variants (column 24, 
# "SNP_ID_CURRENT") per trait (based on EFO ID)
rule count_variantids_disease_specific:
    input: "results/{efo_id}-associations_ontology-annotated.tsv"
    output: "results/{efo_id}_counts_snpidcurrent.tsv"
    shell: "cat {input} | cut -f 24 | sort | uniq -c | sort -k1 -n -r > {output}"

# Count the number of GWAS catalog variants (column 24, 
# "SNP_ID_CURRENT") per trait (based on EFO ID)
rule count_variants_disease_specific:
    input: "results/{efo_id}-associations_ontology-annotated.tsv"
    output: "results/{efo_id}_counts_snps.tsv"
    shell: "cat {input} | cut -f 22 | sort | uniq -c | sort -k1 -n -r > {output}"

# Count the number of GWAS catalog reported genes (column 14, 
# "REPORTED GENE(S)") per trait (based on EFO ID)
rule count_genes_disease_specific:
    input: "results/{efo_id}-associations_ontology-annotated.tsv"
    output: "results/{efo_id}_counts_reportedgenes.tsv"
    shell: "cat {input} | cut -f 14 | sort | uniq -c | sort -k1 -n -r > {output}"

rule get_all_counts:
    input: expand("results/{efo_id}_counts_study.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_snps.tsv",efo_id=EFO_IDS),
           expand("results/{efo_id}_counts_reportedgenes.tsv",efo_id=EFO_IDS)

# We select all studies that have only one mapped trait URI (EFO ID)
# and for which this one EFO ID is within the list of manually obtained 
# Autoimmune disease EFO IDs
rule select_autoimmune_study_by_efo:
    input: "data/autoimmune_efo_ids.tsv",
           "data/gwas_catalog/2022/05/23/gwas-catalog-studies_ontology-annotated.tsv"
    output: "results/autoimmune_gwas-catalog-studies_ontology-annotated.tsv"
    run:
        efo_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                efo_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split("\t")
                 mapped_trait = s[13]
                 if (not ',' in mapped_trait and mapped_trait.split("/")[-1] in efo_ids) or line[:10] == "DATE ADDED":
                    f_out.write(line)

# Get disease-specific lists
rule get_disease_specific_study_lists:
    input: "results/autoimmune_gwas-catalog-studies_ontology-annotated.tsv"
    output: "results/{efo_id}-studies_ontology-annotated.tsv"
    run: 
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split("\t")
                 mapped_trait = s[13]
                 if mapped_trait.split("/")[-1] == wildcards.efo_id or line[:10] == "DATE ADDED":
                    f_out.write(line)

# Count the number of GWAS catalog study accessions (column 37, 
# "STUDY ACCESSION") per trait (based on EFO ID)
rule count_studies_disease_specific_study:
    input: "results/{efo_id}-studies_ontology-annotated.tsv"
    output: "results/{efo_id}_studycounts_study.tsv"
    shell: "cat {input} | cut -f 15 | sort | uniq -c | sort -k1 -n -r > {output}"

rule get_all_counts_study:
    input: expand("results/{efo_id}_studycounts_study.tsv",efo_id=EFO_IDS)

# Extract PGS scores for the EFO lists
rule get_autoimmune_pgs_list:
    input: "data/autoimmune_efo_ids.tsv",
           "data/pgs_catalog/pgs_all_metadata_scores.csv" 
    output: "results/autoimmune_pgs_all_metadata_scores.csv"
    run:
        efo_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                efo_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_efo = s[4]
                 if (not "|" in pgs_efo and pgs_efo in efo_ids) or line[:15] == "Polygenic Score":
                    f_out.write(line)    

# Autoimmune-specific PGS scores for a given EFO ID
# I edited "results/autoimmune_pgs_all_metadata_scores_edited.csv" in Excel to 
# remove "," within the first columns in those cases where content was escaped 
# with "", otherwise the cut with "," doesn't work. 
rule get_disease_specific_autoimmune_pgs_list:
    input: "results/autoimmune_pgs_all_metadata_scores_edited.csv" 
    output: "results/{efo_id}_pgs_all_metadata_scores.csv"
    run: 
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 mapped_trait = s[4]
                 if mapped_trait == wildcards.efo_id or line[:15] == "Polygenic Score":
                    f_out.write(line)   

rule disease_specific_pgs_ids:
    input: "results/{efo_id}_pgs_all_metadata_scores.csv"
    output: "results/{efo_id}_pgsids.tsv"
    shell: "cat {input} | grep -v 'Polygenic Score (PGS) ID' | cut -f 1 -d ',' | sort | uniq  > {output} & true"

rule disease_specific_pgs_pubs:
    input: "results/{efo_id}_pgs_all_metadata_scores.csv"
    output: "results/{efo_id}_pgs_pubcount.tsv"
    shell: "cat {input} | grep -v 'Polygenic Score (PGS) ID' | cut -f 12 -d ',' | sort | uniq  > {output} & true"

rule disease_specific_pgseval_pubs:
    input: "results/{efo_id}_pgs_all_metadata_performance_metrics.csv"
    output: "results/{efo_id}_pgseval_pubcount.tsv"
    shell: "cat {input} | grep -v 'PGS Performance Metric (PPM) ID' | cut -f 4 -d ',' | sort | uniq  > {output} & true"

rule disease_specific_combinedpgs_pubs:
    input: "results/{efo_id}_pgs_all_metadata_scores.csv",
           "results/{efo_id}_pgs_all_metadata_performance_metrics.csv"
    output: "results/{efo_id}_pgscombined_pubcount.tsv"
    shell: "cat {input[0]} | grep -v 'Polygenic Score (PGS) ID' | cut -f 12 -d ',' > {output}.tmp & true; " + \
           "cat {input[1]} | grep -v 'PGS Performance Metric (PPM) ID' | cut -f 4 -d ',' >> {output}.tmp & true; " + \
           "cat {output}.tmp | sort | uniq -c > {output}; " + \
           "rm {output}.tmp; "

rule disease_specific_get_autoimmune_score_development_list:
    input: "results/{efo_id}_pgsids.tsv",
           "results/autoimmune_pgs_all_metadata_score_development_samples.csv"
    output: "results/{efo_id}_pgs_all_metadata_score_development_samples.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[0]
                 if (pgs_id in pgs_ids) or line[:15] == "Polygenic Score":
                    f_out.write(line)  

rule disease_specific_get_autoimmune_score_performance_list:
    input: "results/{efo_id}_pgsids.tsv",
           "results/autoimmune_pgs_all_metadata_performance_metrics.csv" 
    output: "results/{efo_id}_pgs_all_metadata_performance_metrics.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[1]
                 if (pgs_id in pgs_ids) or line[:22] == "PGS Performance Metric":
                    f_out.write(line)  

rule disease_specific_autoimmune_score_evaluation_list:
    input: "results/{efo_id}_pgsids.tsv",
           "results/autoimmune_pgs_all_metadata_evaluation_sample_sets.csv" 
    output: "results/{efo_id}_pgs_all_metadata_evaluation_sample_sets.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[1]
                 if pgs_id in pgs_ids or line[:14] == "PGS Sample Set":
                    f_out.write(line)  

rule count_pgs_evals:
    input: "results/{efo_id}_pgs_all_metadata_performance_metrics.csv"
    output: "results/{efo_id}_pgs_performance_metrics_count.csv"
    shell: "cat {input} | cut -f 1 -d ',' | sort | uniq -c > {output}"

rule count_disease_specific_pgs_studies:
    input: "results/{efo_id}_pgs_all_metadata_scores.csv"
    output: "results/{efo_id}_pgs_study_count.csv"
    shell: "cat {input} | cut -f 12 | sort | uniq -c > {output}"

# Count the number of associations per trait (based on EFO ID)
rule count_pgs_per_disease:
    input: "results/autoimmune_pgs_all_metadata_scores.csv"
    output: "results/autoimmune_pgs_count.tsv"
    shell: "cat {input} | cut -f 5 -d ',' | sort | uniq -c | sort -k1 -n -r > {output}"

# Count the number of associations per trait (based on EFO ID)
rule make_pgs_id_list:
    input: "results/autoimmune_pgs_all_metadata_scores.csv"
    output: "results/autoimmune_pgsids.tsv"
    shell: "cat {input} | cut -f 1 -d ',' | sort | uniq > {output}"

# Extract PGS development samples for autoimmune PGS
rule get_autoimmune_score_development_list:
    input: "results/autoimmune_pgsids.tsv",
           "data/pgs_catalog/pgs_all_metadata_score_development_samples.csv" 
    output: "results/autoimmune_pgs_all_metadata_score_development_samples.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[0]
                 if pgs_id in pgs_ids or line[:15] == "Polygenic Score":
                    f_out.write(line)

# Extract PGS development samples for autoimmune PGS
rule get_autoimmune_score_performance_list:
    input: "results/autoimmune_pgsids.tsv",
           "data/pgs_catalog/pgs_all_metadata_performance_metrics.csv" 
    output: "results/autoimmune_pgs_all_metadata_performance_metrics.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[1]
                 if pgs_id in pgs_ids or line[:15] == "PGS Performance":
                    f_out.write(line)  

# Extract PGS development samples for autoimmune PGS
rule get_autoimmune_score_evaluation_list:
    input: "results/autoimmune_pgsids.tsv",
           "data/pgs_catalog/pgs_all_metadata_evaluation_sample_sets.csv" 
    output: "results/autoimmune_pgs_all_metadata_evaluation_sample_sets.csv"
    run:
        pgs_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                pgs_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.split(",")
                 pgs_id = s[1]
                 if pgs_id in pgs_ids or line[:14] == "PGS Sample Set":
                    f_out.write(line)    

# Extract autoimmune EFO lines from master file
rule get_autoimmune_ukb:
    input: "data/autoimmune_efo_ids.tsv",
           "data/ukb/EFO-UKB-mappings/UK_Biobank_master_file.tsv"
    output: "results/autoimmune_UK_Biobank_master_file.tsv"
    run:
        efo_ids = []
        with open(input[0],"r") as f_in:
            for line in f_in:
                efo_ids.append(line.strip("\n"))
        with open(input[1],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                 s = line.strip("\n").split("\t")
                 assert(len(s)==7)
                 ukb_efo = s[2]
                 if (not "|" in ukb_efo and (ukb_efo in efo_ids) or line[:5] == "ZOOMA"):
                    f_out.write(line) 