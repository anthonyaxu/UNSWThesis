drop view Simple_Mutations_ID;
create or replace view Simple_Mutations_ID
as
    select distinct on (mutation_id)
        id, mutation_id
    from Simple_Somatic_Mutations
;

drop view Mutation_Samples cascade;
create or replace view Mutation_Samples
as
    select mutation_id, sample_id, concat('chr' || chromosome || ':' || chromosome_start || ':' || chromosome_end) as chromosome_seq
    from Simple_Somatic_Mutations
    group by mutation_id, sample_id, chromosome, chromosome_start, chromosome_end
;

create or replace view Mutation_Samples_Count
as
    select mutation_id, chromosome_seq, string_agg(sample_id, ', ' order by sample_id) as sample_id, count(sample_id) as num_samples
    from Mutation_Samples
    group by mutation_id, chromosome_seq
;

drop view Sample_Mutations cascade;
create or replace view Sample_Mutations
as
    select sample_id, mutation_id
    from Simple_Somatic_Mutations
    group by sample_id, mutation_id
;

drop view Sample_Mutations_Count cascade;
create or replace view Sample_Mutations_Count
as
    select sample_id, string_agg(mutation_id, ', ' order by mutation_id) as mutation_id, count(mutation_id) as num_mutations
    from Sample_Mutations
    group by sample_id
    order by num_mutations desc
;

drop view Copy_Number_Sequences cascade;
create or replace view Copy_Number_Sequences
as
    select sample_id, concat('chr' || chromosome || ':' || chromosome_start || ':' || chromosome_end) as seq
    from Copy_Number_Somatic_Mutations
    group by sample_id, seq
    order by sample_id, seq
;

create or replace view Copy_Number_Sequence_Count
as
    select sample_id, string_agg(seq, ', ' order by seq) as seq, count(seq) as num_cnv
    from Copy_Number_Sequences
    group by sample_id
    order by num_cnv desc
;

create or replace view Copy_Number_Unique
as
    select seq
    from Copy_Number_Sequences
    group by seq
    order by seq
;