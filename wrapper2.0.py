
# Importing necessary libraries
import os
import sys
import glob
import argparse
import pandas as pd


# Argument Parsing
def parse_arguments():
    parser = argparse.ArgumentParser(description="Bioinformatics Pipeline Project")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing paired-end FASTQ files (files named *_1.fastq and *_2.fastq)")
    parser.add_argument("--dataset", action='store_true',
                        help="If set, use the sample test data (subset of reads)")
    return parser.parse_args()

# Directory Setup
def create_output_directory():
    out_dir = "PipelineProject_Hillary_Dapsauski"  # Output directory name
    if not os.path.exists(out_dir): 
        os.makedirs(out_dir)
    os.chdir(out_dir)  # All subsequent output will be written here
    return os.getcwd() # Return the full path to the output directory


# Step 1: Download HCMV Genome and CDS Features
def download_and_prepare_cds(log_handle):
    # Check if function is running 
    # log_handle.write("Step 1: Downloading HCMV genome and CDS features...\n")

    # Use NCBI datasets to download the virus genome (including CDS and genome data)
    download_command = "datasets download virus genome accession NC_006273.2 --include cds,genome"
    os.system(download_command)

    # Unzip the downloaded dataset (assumes file is named ncbi_dataset.zip)
    os.system("unzip -o ncbi_dataset.zip")

    # Check if the data is downloaded and unzipped
    # log_handle.write("Downloaded and unzipped HCMV genome data.\n")

# Count the number of CDS features in the CDS FASTA file
def count_cds_features(log_handle):
    # Count the number of CDS features in the CDS FASTA file
    cds_file = "ncbi_dataset/data/cds.fna"
    cds_count = 0
    with open(cds_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                cds_count += 1
    
    log_handle.write(f"The HCMV genome (NC_006273.2) has {cds_count} CDS.\n\n")
    return cds_count


# Step 2 and 3: Build kallisto Index and Run Quantification
def run_kallisto(input, log_handle):
    # Check to see if the function is running
    # log_handle.write("Step 2: Building kallisto index and Step 3: Running kallisto quantification...\n")

    # Build kallisto index using the CDS FASTA file
    cds_file = "ncbi_dataset/data/cds.fna"
    index_command = f"kallisto index -i HCMV_index.idx {cds_file}"

    # Check to see if the index command is running
    # log_handle.write(f"Running command: {index_command}\n")
    os.system(index_command)

    # Create a directory for kallisto output
    os.mkdir("kallisto_output")
    # Get the list of paired-end FASTQ files in the input directory
    samples = sorted(glob.glob(os.path.join(input, "*")))
    
    # Paired-end samples are processed in pairs
    for i in range (0,len(samples),2):
        pair1 = samples[i]
        pair2 = samples[i+1]
        name = pair1.split('/')[-1].split('_')[0]

        # Create a directory for each sample
        os.mkdir(f"kallisto_output/{name}")

        # Run kallisto quantification for each sample
        kallisto_command = f"kallisto quant -i HCMV_index.idx -o kallisto_output/{name} -b 10 -t 2 {pair1} {pair2}"
        
        # Check to see if the kallisto command is running
        #log_handle.write(f"Running kallisto for sample {name}:\n  {kallisto_command}\n")
        os.system(kallisto_command)
    
    # Check to see if the kallisto quantification is completed    
    # log_handle.write("Kallisto quantification completed.\n\n")

# Extract TPM statistics from kallisto output
def extract_tpm_stats(samples):
    #Extract minimum, median, mean, and maximum TPM from abundance.tsv for the samples.
    abundance_file = os.path.join("kallisto_output", samples, "abundance.tsv")

    # Use pandas to read the abundance file
    df = pd.read_csv(abundance_file, sep="\t")
    stats = {
        "min_tpm": df["tpm"].min(),
        "med_tpm": df["tpm"].median(),
        "mean_tpm": df["tpm"].mean(),
        "max_tpm": df["tpm"].max()
    }
    return stats

# Define metadata as a dictionary
metadata = {'SRR5660030': ['Donor 1', '2dpi'],'SRR5660033': ['Donor 1', '6dpi'],'SRR5660044': ['Donor 3', '2dpi'],'SRR5660045': ['Donor 3', '6dpi']}

# Write TPM statistics per sample using metadata info 
def write_tpm_statistics(log_handle):
    # Check to see if the function is running
    # log_handle.write("TPM Statistics per sample:\n")

    header = "\t".join(["sample", "donor", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"])
    log_handle.write(header + "\n")

    # Find Kallisto output directories for each sample
    samples = [os.path.basename(d) for d in glob.glob("kallisto_output/*")]

    # Iterate through each sample, extract TPM stats, and write to log
    for sample in samples:
        donor, condition = metadata[sample]
        stats = extract_tpm_stats(sample) 

        # Write the TPM statistics to the log file
        line = f"{sample}\t{donor}\t{condition}\t{stats['min_tpm']}\t{stats['med_tpm']}\t{stats['mean_tpm']}\t{stats['max_tpm']}"
        log_handle.write(line + "\n")

    log_handle.write("\n")


# Step 4: Differential Expression Analysis with Sleuth
def run_sleuth(log_handle):
    # Check to see if the function is running
    # log_handle.write("Step 4: Running sleuth differential expression analysis...\n")
    
    # Create an input file for the sleuth R script
    sleuth_input = "sleuth_input.txt"

    # Write the sample, condition, and path to the kallisto output directory
    with open(sleuth_input, "w") as f:
        header = "\t".join(["sample", "condition", "path"])
        f.write(header + "\n")
        # For demonstration, we list samples from kallisto_output and assign dummy conditions.
        samples = [os.path.basename(d) for d in glob.glob("kallisto_output/*") if os.path.isdir(d)]

        conditions = ["2dpi", "6dpi"] * (len(samples)//2)

        # Write the sample, condition, and path to the input file  
        for sample, condition in zip(samples, conditions):
            sample_path = os.path.abspath(os.path.join("kallisto_output", sample))
            f.write(f"{sample}\t{condition}\t{sample_path}\n")

    # Run the R script (Sleuth.r must be provided in the repo)
    sleuth_command = "Rscript Sleuth.r"
    # Check to see if the sleuth command is running
    # log_handle.write(f"Running command: {sleuth_command}\n")
    os.system(sleuth_command)

    # Check to see if the sleuth analysis is completed
    # log_handle.write("Sleuth analysis completed. Differential expression results were generated (see sleuth output).\n\n")


# Step 5: Bowtie2 Mapping for Read Filtering
def run_bowtie2_mapping(input, log_handle):
    # Check to see if the function is running
    # log_handle.write("Step 5: Mapping reads with Bowtie2...\n")

    # Create a directory for Bowtie2 output
    os.mkdir("bowtie2_output")

    # Build Bowtie2 index for the HCMV genome 
    genome_fasta = "ncbi_dataset/data/genomic.fna"
    bowtie2_index_prefix = "HCMV_bowtie_index"
    bowtie2_build_cmd = f"bowtie2-build {genome_fasta} {bowtie2_index_prefix}"
    log_handle.write(f"Running command: {bowtie2_build_cmd}\n")
    os.system(bowtie2_build_cmd)
    
    # For each paired-end sample, run Bowtie2 and count read pairs before and after filtering.
    sample_files = sorted(glob.glob(os.path.join(input, "*_1.fastq")))
    for fastq1 in sample_files:
        base = os.path.basename(fastq1).replace("_1.fastq", "")
        fastq2 = os.path.join(input, f"{base}_2.fastq")
        # Count reads before mapping
        count_command1 = f"echo $(( $(wc -l < {fastq1}) ))"
        count_command2 = f"echo $(( $(wc -l < {fastq2}) ))"
        before_count1 = int(os.popen(count_command1).read().strip()) # popen allows you to read the output of the command 
        before_count2 = int(os.popen(count_command2).read().strip())
        before_total = min(before_count1, before_count2)
        
        # Run Bowtie2 mapping; output SAM file will be used to count mapped reads.
        mapped_output = f"bowtie2_output/{base}_bowtie2.sam"
        bowtie2_cmd = f"bowtie2 -x {bowtie2_index_prefix} -1 {fastq1} -2 {fastq2} -S {mapped_output}"
        log_handle.write(f"Running Bowtie2 for sample {base}:\n  {bowtie2_cmd}\n")
        os.system(bowtie2_cmd)
        
        # Count mapped reads using samtools (requires samtools installed)
        mapped_count_cmd = f"samtools view -c -F 4 {mapped_output}"

        # Count reads after mapping
        mapped_count = int(os.popen(mapped_count_cmd).read().strip())
        
        log_handle.write(f"{base} had {before_total} read pairs before Bowtie2 filtering and {mapped_count} read pairs after.\n")
    log_handle.write("\n")

# Step 6: Assemble Reads with SPAdes
def run_spades(log_handle):
    # Make a directory for SPAdes output
    os.mkdir("spades_output")

    # Dictionary of donors and their corresponding 2dpi and 6dpi samples
    donors = {
        'Donor1': {
            '2dpi': ['SRR5660030_1.fastq', 'SRR5660030_2.fastq'],
            '6dpi': ['SRR5660033_1.fastq', 'SRR5660033_2.fastq']
        },
        'Donor3': {
            '2dpi': ['SRR5660044_1.fastq', 'SRR5660044_2.fastq'],
            '6dpi': ['SRR5660045_1.fastq', 'SRR5660045_2.fastq']
        }
    }

    # Process each Bowtie2 SAM file in the bowtie2_output directory
    for sam_file in glob.glob(os.path.join("bowtie2_output", "*.sam")):
        # Extract the sample name from the SAM file name
        sample_name = os.path.basename(sam_file).replace("_bowtie2.sam", "")
        
        # Convert the SAM file to a FASTQ file using samtools
        fastq1 = f"bowtie2_output/{sample_name}_1.fastq"
        fastq2 = f"bowtie2_output/{sample_name}_2.fastq"
        samtools_cmd1 = f"samtools fastq -1 {fastq1} -2 {fastq2} {sam_file}"
        # Check to see if the samtools command is running
        # log_handle.write(f"Running command to convert SAM to FASTQ for {sample_name}:\n  {samtools_cmd1}\n")
        os.system(samtools_cmd1)

    # Ensures SPAdes can find the correct FASTQ files for each donor and condition
    for donor, conditions in donors.items():
        fastq_1_files = os.path.join("bowtie2_output", conditions['2dpi'][0])  
        fastq_1_files_2 = os.path.join("bowtie2_output", conditions['6dpi'][0])  
        fastq_2_files = os.path.join("bowtie2_output", conditions['2dpi'][1])  
        fastq_2_files_2 = os.path.join("bowtie2_output", conditions['6dpi'][1]) 
        
        # Log the SPAdes command being run
        # Check to see if the SPAdes command is running
        # log_handle.write(f"Running SPAdes for {donor} with 2dpi and 6dpi conditions...\n")
        log_handle.write(f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 {fastq_1_files} --pe-2 1 {fastq_2_files} --pe-1 2 {fastq_1_files_2}  --pe-2 2 {fastq_2_files_2} -o spades_output/{donor}\n")

        # Run SPAdes with the combined 2dpi and 6dpi reads for each donor
        spades_command = f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 {fastq_1_files} --pe-2 1 {fastq_2_files} --pe-1 2 {fastq_1_files_2} --pe-2 2 {fastq_2_files_2} -o spades_output/{donor}"
        os.system(spades_command)

        # Check to see if the SPAdes command is completed
        # log_handle.write(f"SPAdes completed for {donor}.\n\n")


# Step 7: BLAST Analysis on Longest Contigs
def run_blast_analysis(log_handle):
    # Check to see if the function is running
    # log_handle.write("Step 7: Running BLAST analysis on longest contigs...\n")

    # For each donor, retrieve the SPAdes contig file
    donors = {
        "Donor1": "spades_output/Donor1/contigs.fasta",
        "Donor3": "spades_output/Donor3/contigs.fasta"
    }
    
    # Iterate through the donors and their corresponding SPAdes contig files
    for donor, contigs_file in donors.items():
        # Check to see if the BLAST analysis is running
        # log_handle.write(f"Processing {donor}...\n")
        
        # Retrieve the longest contig from the SPAdes assembly file
        longest_contig = "" # Initialize the longest contig sequence
        max_len = 0 # Initialize the maximum length of the contig
        header = "" # Initialize the header of the longest contig
        seq = "" # Initialize the sequence of the longest contig
        
        # Open the contigs file and read through each line
        with open(contigs_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    if seq and len(seq) > max_len:
                        longest_contig = header + "\n" + seq
                        max_len = len(seq)
                    header = line.strip()
                    seq = ""
                else:
                    seq += line.strip()
            # Final contig check for the last contig in the file
            if seq and len(seq) > max_len:
                longest_contig = header + "\n" + seq

        # Write the longest contig to a temporary FASTA file
        longest_file = f"{donor}_longest_contig.fasta"
        with open(longest_file, "w") as f_out:
            f_out.write(longest_contig)
        
        # Making blast database 
        downloads_command = "datasets download virus genome taxon Betaherpesvirinae --refseq --include genome"
        os.system(downloads_command)

        # Unzip the downloaded dataset (assumes file is named ncbi_dataset.zip)
        os.system("unzip -o ncbi_dataset.zip")
        
        # Create a local BLAST database using the Betaherpesvirinae genome
        database_command = 'makeblastdb -in ncbi_dataset/data/genomic.fna -dbtype nucl -out Betaherpesvirinae_db'
        os.system(database_command)

        # Run BLAST (blastn) against a local Betaherpesvirinae database 
        blast_output = f"{donor}_blast_results.txt"
        blast_cmd = (
            f"blastn -query {longest_file} -db Betaherpesvirinae_db "
            f"-out {blast_output} -max_target_seqs 10 "
            f"-outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle'"
        )
        # Check to see if the BLAST command is running
        # log_handle.write(f"Running BLAST for {donor}:\n  {blast_cmd}\n")
        os.system(blast_cmd)

        # Write BLAST results with header to the log file
        log_handle.write(f"{donor} BLAST results:\n")
        header_line = "\t".join(["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"])
        log_handle.write(header_line + "\n")
        
        with open(blast_output, "r") as blast_file:
            for line in blast_file:
                log_handle.write(line)
        
        log_handle.write("\n")
    # Check to see if the BLAST analysis is completed
    # log_handle.write("BLAST analysis completed.\n\n")

# Main Function
def main():
    args = parse_arguments()
    cwd = os.getcwd()
    input = os.path.join(cwd, args.input_dir)
    # If the --sample flag is set, instruct the user (or script) to use the sample test data
    if args.dataset:
        print("Using sample test data. (Ensure sample FASTQ files are provided in the input directory.)")
    
    # Create and move into the output directory
    outdir = "PipelineProject_Hillary_Dapsauski"
    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # create output directory and move into it
    os.chdir(outdir)


    # Open the log file where all pipeline output and commands will be written
    with open("PipelineProject.log", "w") as log_handle:
        log_handle.write("Pipeline Project Log\n\n")
        
        # Step 1: Download genome and CDS features; count CDS
        download_and_prepare_cds(log_handle)
        count_cds_features(log_handle)
        
        # Step 2 and 3: Build kallisto index and run quantification
        run_kallisto(input, log_handle) 

        # Write TPM statistics per sample
        write_tpm_statistics(log_handle)
        
        # Step 4: Differential expression analysis with sleuth
        run_sleuth(log_handle)
        
        # Step 5: Map reads with Bowtie2 and log read counts before/after filtering
        run_bowtie2_mapping(input, log_handle)
        
        # Step 6: Assemble reads with SPAdes (one assembly per donor)
        run_spades(log_handle)
        
        # Step 7: Identify the best BLAST hits for the longest contig from each assembly
        run_blast_analysis(log_handle)
    

if __name__ == "__main__":
    main()
