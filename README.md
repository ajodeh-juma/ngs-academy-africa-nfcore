- [Introduction](#introduction)
- [Steps](#steps)


## Introduction
This lesson is an introduction to the workflow manager Nextflow, and nf-core, a community effort to collect a curated set of analysis pipelines built using Nextflow.

Nextflow enables scalable and reproducible scientific workflows using software containers such as Docker and Singularity. It allows the adaptation of pipelines written in the most common scripting languages such as R and Python. Nextflow is a Domain Specific Language (DSL) that simplifies the implementation and the deployment of complex parallel and reactive workflows on clouds and clusters.

This lesson motivates the use of Nextflow and nf-core as a development tool for building and sharing computational pipelines that facilitate reproducible (data) science workflows.


## Steps
#### ***Log into the HPC***
From the terminal of your local computer, you can log into the HPC using the
following command line, followed by pressing <ENTER>. You will be prompted to
type in your password. On a Linux system, you can use the `Ctrl-Alt-T` keyboard
shortcut to open a terminal.
```
ssh <train11>@172.16.13.171
```

If using PuTTY, type `172.16.13.171` in the `Host Name (or IP address)` field.


In your `home` directory, follow the steps:

1. Clone the repo 
