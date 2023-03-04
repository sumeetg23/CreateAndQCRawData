# Snakemake workflow for demultiplex and QC of NGS runs

Snakemake workflow to generate demultiplexed FASTQ files, FASTQC and FASTQSCREEN reports for each fastq file, and summarize FASTQ and FASTQSCREEN reports on a per lane basis.

```
Command:
snakemake -s rules/bcl2fastq.smk --cluster-config clusterTime.json --cluster "sbatch --mem={cluster.mem} --cpus-per-task={cluster.cpus_per_task} --account={cluster.account} --partition={cluster.slurm_partition} --output %j.out --wrap" --jobs 100 --latency-wait 60
```
