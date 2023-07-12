# MSK Nextflow Training

## Install Nextflow
helpful link: 
https://www.nextflow.io/index.html

## Prerequisite

- Docker:
  If not installed you can install here: https://www.docker.com/

- Java: 
```
java -version
```

If you do no have Java installed you can install with: 

```
sdk install java 17.0.6-amzn
```

## Install Nextflow
```
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

There may be some permission issues. Please reachout to me if you experience these.

resources: https://www.nextflow.io/ 


## Training Repo

clone the training repo with: 
```
git clone https://github.com/mskcc-omics-workflows/nextflow-training.git
```

## Run Example Workflow


Simple: 
```
nextflow workflow.nf 
```

With Parameter:
```
nextflow workflow.nf --in /path/to/nextflow-training/data/transcriptome.fa
```

## Using Containers and Configs

First try the following: 

```
nextflow container.nf --transcript /Users/ebuehler/Documents/GitHub/nextflow-training/transcriptome.fa 
```

This command should fail, which makes sense as the `salmon` command is not necessarily on our local system. 

So far all of our commands have been avaliable on a local machine, but what if we want to call a command that exists in a docker or singularity image? Often times we'd prefer to use a containerized environment.

Nextflow makes this easy. Instead simply try: 

```
nextflow container.nf --transcript /Users/ebuehler/Documents/GitHub/nextflow-training/transcriptome.fa  -with-docker
```

But how does Nextflow know how to do this does this? 

It turns out that this container path is already configured for us in the `nextflow.config` file.

Additionally, we can forgo the `-with-docker` option by un-commenting `docker.enabled = true` in the `nextflow.config` file. 

Then running works: 

```
nextflow container.nf --transcript /Users/ebuehler/Documents/GitHub/nextflow-training/transcriptome.fa 
```

However, we maybe want to use configs different from this one. 

We can specify a different `nextflow.config` file by doing the following:

```
nextflow container.nf --transcript /Users/ebuehler/Documents/GitHub/nextflow-training/transcriptome.fa  -with-docker -c /Users/ebuehler/Documents/GitHub/nextflow-training/alternative.config
```

On the surface, it appears this command did the same as the last, but it in fact they used different containers. You can see this by comparing the `process.container` argument in `nextflow.config` vs the `alternative.config`. 

## Multiple Processes and Multiple Inputs

In any Bioinformatics pipeline, we'd like to combine these processses together. We can see this by running the following DSL-2 example: 

```
nextflow processes.nf
```

* Note: that Nextflow processes define the execution of asynchronous tasks i.e. they are not executed one after another as they are written in the pipeline script as it would happen in a common imperative programming language.

In other words, the only thing keeping index from running before quantification is the fact that index output is used for quantification. 

We can see how this might work otherwise by commenting out: 

```
quantification(index.out.index, read_pairs_ch)
```

and un-commenting: 
```
// index_dir = "$baseDir/nextflow-training/data/index"
// quantification(index_dir, read_pairs_ch)
```
in `processes.nf` and then re-running the script with: 

```
nextflow processes.nf
```

you will notice both processes start simultaneously. 

- Multiple inputs
In `processes.nf`, if we comment out:

```
params.reads = "$baseDir/data/gut_{1,2}.fq"
```

and un-comment: 

```
// params.reads = "$baseDir/data/*_{1,2}.fq"
```

We can see what it's like to pass multiple inputs into a process. Now instead of just passing one pair of fq files, we pass 2 sets, gut and liver. 

Now quantification reports processing both of these sets: `2 of 2 ✔`. 

Also note that when multiple inputs are passed to a process as Channels they are processed in parrallel: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel

If you uncomment `read_pairs_ch.view()` you can take a look at what the `read_pairs_ch` channel looks like: 

```
[liver, [/Users/ebuehler/Documents/GitHub/nextflow-training/data/liver_1.fq, /Users/ebuehler/Documents/GitHub/nextflow-training/data/liver_2.fq]]
[gut, [/Users/ebuehler/Documents/GitHub/nextflow-training/data/gut_1.fq, /Users/ebuehler/Documents/GitHub/nextflow-training/data/gut_2.fq]]
```
Notice how each pair of files also include a tag that distinguishes. This dictionary type structure is very typical of nextflow Channels, and includes a variety of methods: https://www.nextflow.io/docs/latest/operator.html.

Keep this distinct from the idea of asynchronous processsing of processes.

## Pulling Data from Other Sources
Thus far all of our data has come from this repositories local data directory. 

However, we can use the power of nextflow to pull from online sources as well. 

To do this, try runnning: 
```
nextflow processes_data.nf -c alternative_data.config
```
This script specifies yet another config that tells our script where some github hosted data is. We an then pull that data into our program and combine it with our local data.

Notice how the output now identifies three runs of quantification. This features can be particularly powerful for pulling down test data or data into a variety of environments. 

 
- introduce alignment in regular nextflow

## Alignment, a CCI Example

Now it's time to try an example a little more specific to the needs of CCI. 

To start first install some additional test data by running: 
```
source install_data.sh
```

This will install the folder test_nucleo to our current directory. 

Next try running the following nextflow script:
```
nextflow run alignment.nf -entry test_alignment -c advanced.config
```

There are a few things to note here first we are using a different config file, `advanced.config`. Additionally, there is another parameter named `-entry`, which specifies which workflow within the the nextflow script we'd like to enter. 

We can see these different workflows by taking a closer look at the `alignment.nf` script.

* Note that within this script we notice a few new things. 

First are the include clauses: 
```
include { BWA_MEM      } from './modules/bwa_mem.nf'
include { PICARD_ADDORREPLACEREADGROUPS      } from './modules/picard_addorreplacereadgroups.nf'
```

These allow us to import nextflow processes for use in the workflow defined in `alignment.nf`. This proves useful as it is easy to see how a single nextflow script might get crowded without this functionality.

This script also does more to show off the power of Channel Operators. 


This proves to be powerful example of the uses of nextflow. And from here we have enought to write a bioninformatics pipeline. However, there are still a few weaknesses. 

For one, advanced.config is getting quite bloated. It is true we can specify multiple configs. However, it is not exactly clear what would be the best way to do this. In a similar note, it is also not clear how to best organize our modules and workflows. Or how testing should be handled in a formalized way? 

It seems we are missing some best practice standards and accompanying tools to help us navigate these obstacles. 

Of course, this is something anyone could be solved by an individual team applying good software enginneering practices. Alternatively, it turns out there is already an effort to do this via nf-core.

# Modularizing Alignment in Nf-Core


- lucid chart of nf-core / explain our fork
- pull into pipeline 
To get started, we can install the nf-core python package with the following: 

```
pip install nf-core==2.7.1
```

* Note please install 2.7.1 as this is what materials for this training have been made using. 2.9 is a recent release and I haven't vetted for compatibility.

Creating a pipeline: 
```
git remote add origin https://github.com/mskcc-omics-workflows/pipeline_example.git
git branch -M main 
git push -u origin main
```

Installing a module: 

```
nf-core modules install bwa_mem
nf-core modules install picard_addorreplacereadgroups
nf-core subworkflows --git-remote https://github.com/mskcc-omics-workflows/mskcc-modules.git --branch feature/extractumi install extractumi
```

Other topics of interest:
- Singularity? 
- Running on Cluster? 
- Execution from online
- How to contribute to mskcc-omics-workflows/modules 



Pieces adapted from Sequera: https://github.com/seqeralabs 
