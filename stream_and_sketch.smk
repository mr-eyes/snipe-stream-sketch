import os
import pandas as pd

# Load configuration
configfile: "config.yaml"

# Helper functions to filter and extract relevant information from the dataframe
def load_filtered_df():
    """Dynamically load and filter the dataframe."""
    df = pd.read_csv(config["raw_data_file"], sep="\t", low_memory=False)
    df = df[df["library_strategy"] == 'WGS']
    df['PE1'] = df['fastq_ftp'].apply(lambda x: x.split(';')[0] if isinstance(x, str) else None)
    df['PE2'] = df['fastq_ftp'].apply(lambda x: x.split(';')[1] if isinstance(x, str) and ';' in x else None)
    df['SE'] = df['fastq_ftp'].apply(lambda x: isinstance(x, str) and ';' not in x)
    df = df[df['SE'] == True].head(1)
    return df

def get_pe1(run_id):
    df = load_filtered_df()
    return df[df['run_accession'] == run_id]['PE1'].values[0]

def get_pe2(run_id):
    df = load_filtered_df()
    return df[df['run_accession'] == run_id]['PE2'].values[0]

def is_single_end(run_id):
    df = load_filtered_df()
    return df[df['run_accession'] == run_id]['SE'].values[0]

def get_run_ids():
    df = load_filtered_df()
    return df['run_accession'].tolist()

# Final target rule to ensure all signatures are created
rule all:
    input:
        expand(os.path.join(config["signatures_dir"], "{run_id}.sig"), run_id=get_run_ids())

# Rule for downloading and sketching with sourmash
rule download_and_sketch:
    output:
        sig=os.path.join(config["signatures_dir"], "{run_id}.sig"),
        stdout_log=os.path.join(config["logs_dir"], "{run_id}.out"),
        stderr_log=os.path.join(config["logs_dir"], "{run_id}.err"),
        sketched_seq_log=os.path.join(config["log_number_of_sketched_seqs_dir"], "{run_id}_sketched_seqs.log")
    params:
        pe1=lambda wildcards: get_pe1(wildcards.run_id),
        pe2=lambda wildcards: get_pe2(wildcards.run_id),
        se=lambda wildcards: is_single_end(wildcards.run_id)
    shell:
        """
        if [[ "{params.se}" == "True" ]]; then
            # Single-end case: stream one FASTA file
            curl -s {params.pe1} | zcat | sourmash sketch dna - -p k=51,scaled=10000,abund --name {wildcards.run_id} -o {output.sig} \
            > {output.stdout_log} 2> {output.stderr_log}
        else
            # Paired-end case: concatenate R1 and R2 to simulate a single FASTA stream
            cat <(curl -s {params.pe1} | zcat) <(curl -s {params.pe2} | zcat) | \
            sourmash sketch dna - -p k=51,scaled=10000,abund --name {wildcards.run_id} -o {output.sig} \
            > {output.stdout_log} 2> {output.stderr_log}
        fi

        # Extract the line containing "sequences taken from 1 files"
        grep "sequences taken from 1 files" {output.stderr_log} > {output.sketched_seq_log}
        """
