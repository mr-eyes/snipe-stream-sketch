import pandas as pd

# Load the list of run IDs from the CSV file
run_ids = pd.read_csv("mice_sra_runinfo.csv.gz", sep=',', usecols=["Run"], compression='gzip')["Run"].tolist()
print(f"Number of runs: {len(run_ids)}")


rule all:
    input:
        expand("PE_SIGS_1K/{sample}.sig", sample=run_ids),
        expand("PE_SIGS_10K/{sample}.sig", sample=run_ids),
        expand("CLEANUP_DONE/cleanup_complete_{sample}.txt", sample=run_ids)


rule sra:
    output:
         "sra/{sample}/{sample}"
    threads: 1
    priority: 0
    resources:
        mem_mb=lambda wildcards, attempt: 5 * 1024 * attempt,
        time=lambda wildcards, attempt: 48 * 60 * attempt,
        runtime=lambda wildcards, attempt: 48 * 60 * attempt,
        partition="bmm"
    shell:
        "aws s3 sync --no-sign-request s3://sra-pub-run-odp/sra/{wildcards.sample} sra/{wildcards.sample}/"


rule skt1k:
    input:
        "sra/{sample}/{sample}"
    output:
        "PE_SIGS_1K/{sample}.sig"
    threads: 4
    priority: 10
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        time=lambda wildcards, attempt: 48 * 60 * attempt,
        runtime=lambda wildcards, attempt: 48 * 60 * attempt,
        partition="bmm"
    params:
        sample_name="{sample}"
    shell:
        """
        fasterq-dump --skip-technical --fasta-unsorted --threads {threads} --bufsize 1000MB --curcache 10000MB --mem {resources.mem_mb} {input} --stdout | sourmash sketch dna - -p k=51,scaled=1000,abund -o {output} --name {params.sample_name}
        """


rule sk10k:
    input:
        sig="PE_SIGS_1K/{sample}.sig"
    output:
        "PE_SIGS_10K/{sample}.sig"
    params:
        sample_name="{sample}"
    priority: 40
    resources:
        mem_mb=lambda wildcards, attempt: 5 * 1024 * attempt,
        time=lambda wildcards, attempt: 4 * 60 * attempt,
        runtime=lambda wildcards, attempt: 4 * 60 * attempt,
        partition="low2"
    shell:
        "sourmash signature downsample -q {input.sig} --scaled 10000 -k 51 -o {output}"


rule clean:
    priority: 100
    input:
        sig1="PE_SIGS_1K/{sample}.sig",
        sig2="PE_SIGS_10K/{sample}.sig"
    output:
        "CLEANUP_DONE/cleanup_complete_{sample}.txt"
    shell:
        """
        rm -rf sra/{wildcards.sample}/
        touch {output}
        """
