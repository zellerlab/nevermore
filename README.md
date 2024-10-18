## ZellerLab Nevermore

This is forked from Christians Nevermore Pipeline at EMBL.  
[https://github.com/cschu/nevermore](https://github.com/cschu/nevermore)  
This Pipeline is currently still a work in progress 
and will be slimmed down to accomodate the Groups needs.  
As such there are still many bits which can be slimmed down to suit our specific needs.
We have plenty of ideas on how to improve this too.  
  
In the mean time, create Github issues to complain!

### What does this pipeline do ? - The (very) high level overview
It processes Shotgun Metagenomic files. This is a basic general workflow.  
Data Preprocessing:  
- input file handling
- Quality Control
- Human Decontamination
- Stats collection

Read Alignment to the GMGC:  
- either alignment is wrapped by Gffquant OR
- alignment and gffquant are performed sequentially

GFFQuant:  
- Functional Annotations from the GMGC are transferred to aligned Shotgun Reads

Running MoTus3, humann3 or MetaPhlan4 Taxonomic Profiling  

### How do I get this pipeline?
Currenlty I recommend cloning your own local version of this repo.
This gives you the updated config and param files and 
makes version control a bit more stable while the repo is a WIP.
To do so (for example on the SHARK cluster at LUMC), run:  

```
cd /exports/archive/lucid-grpzeller-primary/$USER
git clone -b test1 https://github.com/zellerlab/nevermore.git
cd nevermore
```
You can find the config files under `nevermore/config/`  

### How do I run this?
Edit either the `run_lumc.sh` or `run_embl.sh` file to suit your needs.  
- Give your run a sensible job_name (for slurm)
- pick suitable slurm runtimes for the main slurm job
- Adapt your WORKFLOW path
- Set a suitable WORKDIR path
- Set suitable params and config files
- adapt the queue keyword in the `/nevermore/config/run.config` file to the cluster

Edit the `params.yml` and `config.yml` files to suit your needs  
  
IMPORTANT: After your run has successfully completed:  
  
your files will the in the output_dir specified in your params file.  
Copy it to your archived parition, resolving symlinks using:  
`cp -Lr outdir archived_distination_dir`  

save a copy of params.yml and config.yml together with your run!  
`cp {params.yml,run.config} archived destination_dir` 

Submit the job with:  
`sbatch run_lumc.sh`  

### Requirements
This pipeline is meant to be run on the cluster environments at either EMBL or LUMC.
These requriements are satisfied on both clusters and can be handled in your run.sh script

Formal requirements include:  
- nextflow>=22.10.6 (it is written for DSL2)
- python>=3.7
- gcc (gcc.8.5.0 works)
- internet access to pull containers

## Happy Pipelining!