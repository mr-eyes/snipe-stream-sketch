import pandas as pd
configfile: "config.yaml"


# Load and filter dataframe, then save relevant columns for later use
df = pd.read_csv(config["raw_data_file"], sep="\t", low_memory=False)
df = df[df["library_strategy"].isin(['WGS', 'WXS'])]
df = df[df["instrument_platform"].isin(['ILLUMINA'])]
df = df.dropna(subset=['fastq_ftp'])

# Split PE1, PE2, SE
df['PE1'] = df['fastq_ftp'].apply(lambda x: x.split(';')[0] if isinstance(x, str) else None)
df['PE2'] = df['fastq_ftp'].apply(lambda x: x.split(';')[1] if isinstance(x, str) and ';' in x else None)
df['SE'] = df['fastq_ftp'].apply(lambda x: isinstance(x, str) and ';' not in x)

# Save processed data
df[['run_accession', 'PE1', 'PE2', 'SE']].to_csv("preprocessed_runs.csv", index=False)

runs_df = df.copy()

run_ids = runs_df['run_accession'].tolist()

print(f"Number of runs: {len(runs_df)}")    



# Final target rule to ensure all signatures are created
rule all:
    input:
        expand(os.path.join(config["signatures_dir"], "{run_id}.sig"), run_id=run_ids)

rule download_and_sketch:
    output:
        sig=os.path.join(config["signatures_dir"], "{run_id}.sig"),
        stdout_log=os.path.join(config["logs_dir"], "{run_id}.out"),
        stderr_log=os.path.join(config["logs_dir"], "{run_id}.err"),
        sketched_seq_log=os.path.join(config["log_number_of_sketched_seqs_dir"], "{run_id}_sketched_seqs.log")
    params:
        pe1=lambda wildcards: runs_df[runs_df['run_accession'] == wildcards.run_id]['PE1'].values[0],
        pe2=lambda wildcards: runs_df[runs_df['run_accession'] == wildcards.run_id]['PE2'].values[0],
        se=lambda wildcards: runs_df[runs_df['run_accession'] == wildcards.run_id]['SE'].values[0]
    shell:
        """
        TMPDIR=$(mktemp -d)

        if [[ "{params.se}" == "True" ]]; then
            # Single-end case: download and process
            curl -L -C - --retry 10 --retry-connrefused --retry-delay 5 --max-time 600 --limit-rate 1M \
                {params.pe1} --output $TMPDIR/pe1.fastq.gz
            gunzip -c $TMPDIR/pe1.fastq.gz | sourmash sketch dna - -p k=51,scaled=10000,abund --name {wildcards.run_id} -o {output.sig} \
            > {output.stdout_log} 2> {output.stderr_log}
        else
            # Paired-end case: download both files and process
            curl -L -C - --retry 10 --retry-connrefused --retry-delay 5 --max-time 600 --limit-rate 1M \
                {params.pe1} --output $TMPDIR/pe1.fastq.gz
            curl -L -C - --retry 10 --retry-connrefused --retry-delay 5 --max-time 600 --limit-rate 1M \
                {params.pe2} --output $TMPDIR/pe2.fastq.gz
            gunzip -c $TMPDIR/pe1.fastq.gz $TMPDIR/pe2.fastq.gz | sourmash sketch dna - -p k=51,scaled=10000,abund --name {wildcards.run_id} -o {output.sig} \
            > {output.stdout_log} 2> {output.stderr_log}
        fi

        # Extract the line containing "sequences taken from 1 files"
        grep "sequences taken from 1 files" {output.stderr_log} > {output.sketched_seq_log}

        # Clean up temporary files
        rm -rf $TMPDIR
        """
