- [Introduction](#introduction)
- [Login](#login)
- [Prep](#prep)
- [Your First script](#your-first-script)
- [Simple RNA-Seq pipeline](#simple-rna-seq-pipeline)



## Introduction
This lesson is an introduction to the workflow manager Nextflow, and nf-core, a community effort to collect a curated set of analysis pipelines built using Nextflow.

Nextflow enables scalable and reproducible scientific workflows using software containers such as Docker and Singularity. It allows the adaptation of pipelines written in the most common scripting languages such as R and Python. Nextflow is a Domain Specific Language (DSL) that simplifies the implementation and the deployment of complex parallel and reactive workflows on clouds and clusters.

This lesson motivates the use of Nextflow and nf-core as a development tool for building and sharing computational pipelines that facilitate reproducible (data) science workflows.


## Login
From the terminal of your local computer, you can log into the HPC using the
following command line, followed by pressing <ENTER>. You will be prompted to
type in your password. On a Linux system, you can use the `Ctrl-Alt-T` keyboard
shortcut to open a terminal.
```
ssh <train11>@172.16.13.171
```

If using PuTTY, type `172.16.13.171` in the `Host Name (or IP address)` field
and `open` the program. Login with user name when prompted and key in your password.

## Prep
In your `home` directory, follow the steps:

Clone the repo in your `home` directory
    ```
    git clone https://github.com/ajodeh-juma/ngs-academy-africa-nfcore.git
    ```
## Your first script
1. Open your first nextflow script `wc.nf` using your favourite text editor
   (`nano` or `vim`)
2. Run the script using `nextflow`
    ```
    nextflow run wc.nf
    ```
2. Create a `process` in the script to `print` the number of reads in the input
   file provided. Ensure that you capture the `output` as `stdout`

    **Quiz:** *How many reads are in the input file?*

    ---
    <details close>
    <summary>Answer</summary>
    
    </details>

    ---

## Simple RNA-Seq pipeline

### Using Conda/Bioconda

1. Create a `conda` environment
    ```
    conda env create -f environment.yml
    ```
    The `environment.yml` file has all the required tools/software and
    dependencies for the simple pipleine that we will run. (You can have a
    preview of the file)

2. Activate the conda environment
    ```
    conda activate rnaseq-env
    ```

3. Run the script using `conda` profile
    ```
    nextflow run main.nf -profile conda
    ```

4. Deactivate the `conda` environment
    ```
    conda deactivate
    ```

### Using Docker

1. Build a Docker image
    ```
    docker build -t rnaseq-image .
    ```
    This may take a couple of minutes

2. Test container by looking at the `Salmon` version
    ```
    docker run rnaseq-image salmon --version
    ```

2. Run the script using `docker`
    ```
    nextflow run main.nf -with-docker rnaseq-image
    ```


## Challenge
1. In your `workflow`:
    
    (a). Add a `process` that preprocesses the raw reads using `fastp` and use
    the preprocessed reads as input for the `quantification` step with `salmon`

    (b). Add a `process` that counts the number of preprocessed reads. Print the
    output in `stdout`.
    
2. Run the `workflow` using the `conda` profile.


