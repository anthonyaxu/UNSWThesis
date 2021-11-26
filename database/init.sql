-- sudo service postgresql start
-- sudo -u postgres psql
-- psql -U anthony thesis
-- psql -d thesis -U anthony -f [FILE]

drop table if exists Simple_Somatic_Mutations_Raw CASCADE;
create table Simple_Somatic_Mutations_Raw (
    icgc_mutation_id                text,
    icgc_donor_id                   text,
    project_code                    text,
    icgc_specimen_id                text,
    icgc_sample_id                  text,
    matched_icgc_sample_id          text,
    submitted_sample_id             text,
    submitted_matched_sample_id     text,
    chromosome                      text,
    chromosome_start                integer,
    chromosome_end                  integer,
    chromosome_strand               integer,
    assembly_version                text,
    mutation_type                   text,
    reference_genome_allele         text,
    mutated_from_allele             text,
    mutated_to_allele               text,
    quality_score                   float,
    probability                     float,   
    total_read_count                integer,
    mutant_allele_read_count        integer,
    verification_status             text,
    verification_platform           text,
    biological_validation_status    text,
    biological_validation_platform  text,
    consequence_type                text,    
    aa_mutation                     text,
    cds_mutation                    text,
    gene_affected                   text,
    transcript_affected             text,
    gene_build_version              text,
    platform                        text,
    experimental_protocol           text,
    sequencing_strategy             text,
    base_calling_algorithm          text,
    alignment_algorithm             text,
    variation_calling_algorithm     text,
    other_analysis_algorithm        text,
    seq_coverage                    text,
    raw_data_repository             text,
    raw_data_accession              text,
    initial_data_release_date       text
);
\copy Simple_Somatic_Mutations_Raw from 'simple_somatic_mutation.open.tsv' delimiter E'\t' csv header;
alter table Simple_Somatic_Mutations_Raw add column id SERIAL PRIMARY KEY;
drop table if exists Simple_Somatic_Mutations CASCADE;
create table Simple_Somatic_Mutations
as
    select id,
           icgc_mutation_id as mutation_id,
           icgc_specimen_id as specimen_id,
           icgc_sample_id as sample_id,
           chromosome,
           chromosome_start,
           chromosome_end,
           mutation_type,
           mutated_from_allele,
           mutated_to_allele,
           total_read_count,
           mutant_allele_read_count,
           consequence_type,
           gene_affected,
           transcript_affected
    from Simple_Somatic_Mutations_Raw
    order by id
;
\copy Simple_Somatic_Mutations to 'simple_somatic_mutations.csv' csv header;


drop table if exists Structural_Somatic_Mutations_Raw CASCADE;
create table Structural_Somatic_Mutations_Raw (
    icgc_donor_id                       text,
    project_code                        text,
    icgc_specimen_id                    text,
    icgc_sample_id                      text,
    submitted_sample_id                 text,
    submitted_matched_sample_id         text,
    variant_type                        text,
    sv_id                               text,
    placement                           text,
    annotation                          text,
    interpreted_annotation              text,
    chr_from                            text,
    chr_from_bkpt                       float,   
    chr_from_strand                     integer,
    chr_from_range                      text,
    chr_from_flanking_seq               text,
    chr_to                              text,
    chr_to_bkpt                         float,   
    chr_to_strand                       integer,
    chr_to_range                        text,
    chr_to_flanking_seq                 text,
    assembly_version                    text,
    sequencing_strategy                 text,
    microhomology_sequence              text,
    non_templated_sequence              text,
    evidence                            text,
    quality_score                       float,
    probability                         float,
    zygosity                            text,
    verification_status                 text,
    verification_platform               text,
    gene_affected_by_bkpt_from          text,
    gene_affected_by_bkpt_to	        text,
    transcript_affected_by_bkpt_from    text,
    transcript_affected_by_bkpt_to	    text,
    bkpt_from_context	                text,
    bkpt_to_context	                    text,
    gene_build_version	                text,
    platform	                        text,
    experimental_protocol	            text,
    base_calling_algorithm	            text,
    alignment_algorithm	                text,
    variation_calling_algorithm	        text,
    other_analysis_algorithm	        text,
    seq_coverage	                    float,
    raw_data_repository	                text,
    raw_data_accession                  text
);
\copy Structural_Somatic_Mutations_Raw from 'structural_somatic_mutation.tsv' delimiter E'\t' csv header;
alter table Structural_Somatic_Mutations_Raw add column id SERIAL PRIMARY KEY;
drop table if exists Structural_Somatic_Mutations CASCADE;
create table Structural_Somatic_Mutations
as
    select id,
           icgc_donor_id as donor_id,
           icgc_specimen_id as specimen_id,
           icgc_sample_id as sample_id,
           submitted_sample_id,
           submitted_matched_sample_id,
           variant_type,
           sv_id,
           annotation,
           interpreted_annotation,
           chr_from,
           chr_from_bkpt,
           chr_to,
           chr_to_bkpt,
           microhomology_sequence,
           non_templated_sequence,
           evidence,
           seq_coverage
    from Structural_Somatic_Mutations_Raw
    order by id
;
\copy Structural_Somatic_Mutations to 'structural_somatic_mutations.csv' csv header;


drop table if exists Copy_Number_Somatic_Mutations_Raw CASCADE;
create table Copy_Number_Somatic_Mutations_Raw (
    icgc_donor_id	                text,
    project_code                    text,
    icgc_specimen_id	            text,
    icgc_sample_id	                text,
    matched_icgc_sample_id	        text,
    submitted_sample_id	            text,
    submitted_matched_sample_id	    text,
    mutation_type	                text,
    copy_number	                    float,
    segment_mean	                float,
    segment_median	                float,
    chromosome	                    text,
    chromosome_start	            integer,
    chromosome_end	                integer,
    assembly_version	            text,
    chromosome_start_range	        text,
    chromosome_end_range	        text,
    start_probe_id	                integer,
    end_probe_id	                integer,
    sequencing_strategy	            text,
    quality_score	                integer,
    probability	                    float,
    is_annotated	                text,
    verification_status	            text,
    verification_platform	        text,
    gene_affected	                text,
    transcript_affected	            text,
    gene_build_version	            text,
    platform	                    text,
    experimental_protocol	        text,
    base_calling_algorithm	        text,
    alignment_algorithm	            text,
    variation_calling_algorithm	    text,
    other_analysis_algorithm	    text,
    seq_coverage	                float,
    raw_data_repository	            text,
    raw_data_accession              text
);
\copy Copy_Number_Somatic_Mutations_Raw from 'copy_number_somatic_mutation.tsv' delimiter E'\t' csv header;
alter table Copy_Number_Somatic_Mutations_Raw add column id SERIAL PRIMARY KEY;
drop table if exists Copy_Number_Somatic_Mutations CASCADE;
create table Copy_Number_Somatic_Mutations
as
    select id,
           icgc_donor_id as donor_id,
           icgc_specimen_id as specimen_id,
           icgc_sample_id as sample_id,
           matched_icgc_sample_id as matched_sample_id,
           submitted_sample_id,
           submitted_matched_sample_id,
           mutation_type,
           copy_number,
           segment_mean,
           segment_median,
           chromosome,
           chromosome_start,
           chromosome_end,
           chromosome_start_range,
           chromosome_end_range,
           seq_coverage
    from Copy_Number_Somatic_Mutations_Raw
    order by id
;

\copy Copy_Number_Somatic_Mutations to 'copy_number_somatic_mutations.csv' csv header;