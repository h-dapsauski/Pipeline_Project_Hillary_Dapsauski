# Pipeline_Project_Hillary_Dapsauski
# Bioinformatics Pipeline Documentation

## Required Dependencies

To utilize and run the pipeline, install the following dependencies:

- **Python** - [Download](https://www.python.org/downloads/)
- **Biopython** - [Download](https://biopython.org/wiki/Download)
- **Kallisto** - [Installation Guide](https://github.com/pachterlab/kallisto/blob/master/INSTALL.md)
- **Sleuth** - [About Sleuth](https://pachterlab.github.io/sleuth/about)
- **Bowtie2** - [Installation Guide](https://www.metagenomics.wiki/tools/bowtie2/install)
- **SPAdes** - [GitHub](https://github.com/ablab/spades)
- **BLAST+** - [NCBI Guide](https://www.ncbi.nlm.nih.gov/books/NBK569861/)

---

## Step 1: Downloading and Preparing Data

### Part 1: Downloading SRA Files

The necessary SRA files were downloaded using the `wget` command in the terminal.

1. Obtain the SRA Run ID from the provided web links.
2. Navigate to the NCBI SRA pages and select the **Run ID**.
3. Click the **Data access** tab to retrieve the AWS SRA link.
4. Use `wget` to download the data:

```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```

### Part 2: Converting SRA to FASTQ

The downloaded SRA files were converted into paired-end FASTQ files using `fasterq-dump`:

```bash
fasterq-dump SRR5660030
fasterq-dump SRR5660033
fasterq-dump SRR5660044
fasterq-dump SRR5660045
```

### Part 3: Creating a Sample Dataset

To create a smaller dataset for testing, a new directory **sample_data** was created. The first 40,000 lines of each FASTQ file were extracted and stored in this directory:

```bash
mkdir sample_data && for file in SRR*.fastq; do head -n 40000 "$file" > "sample_data/$file" & done; wait
```

### Part 4: Organizing the Original Data

All original FASTQ files were moved to a separate directory named **original_data**:

```bash
mkdir original_data && mv SRR*.fastq original_data/
```

---

This documentation serves as an initial guide for setting up and preparing data for the bioinformatics pipeline. Additional steps will follow for transcriptome analysis, mapping, assembly, and annotation.

