# MSK Nextflow Training

## Set-up
helpful link: 
https://www.nextflow.io/index.html

### Prerequisites

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

### Install Nextflow
```
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

There may be some permission issues. Please reachout to me if you experience these.

resources: https://www.nextflow.io/ 


### Training Repo

clone the training repo with: 
```
git clone https://github.com/mskcc-omics-workflows/nextflow-training.git
```

## Getting Started with Nextflow

### The Basics 
Now that everything is setup, we can get started by running a simple nextflow script: 
```
nextflow workflow.nf 
```

With Parameter:
```
nextflow workflow.nf --in $(PWD)/data/transcriptome.fa
```

### Using Containers and Configs

First try the following: 

```
nextflow container.nf --transcript $(PWD)/data/transcriptome.fa 
```

This command should fail, which makes sense as the `salmon` command is not necessarily on our local system. If it succeeded, it means you happened to have salmon already installed on your local system.

So far all of our commands have been avaliable on a local machine, but what if we want to call a command that exists in a docker or singularity image? Often times we'd prefer to use a containerized environment.

Nextflow makes this easy. Instead simply try: 

```
nextflow container.nf --transcript $(PWD)/data/transcriptome.fa  -with-docker
```

But how does Nextflow know how to do this does this? 

It turns out that this container path is already configured for us in the `nextflow.config` file.

Additionally, we can forgo the `-with-docker` option and specify a variety of docker configurations. For example try running: 

```
nextflow container.nf --transcript $(PWD)/data/transcriptome.fa -c alternative_docker.config
```
`alternative_docker.config` sets the `docker.enabled` parameter so we don't have to add the `-with-docker` as well as some other docker parameters. 


Note that although we specified a different config file than `nextflow.config`, Nextflow still finds `nextflow.config` and combines it's parameters with the specified configs.

It does this on a directory level so if we'd like to ignore `nextflow.config`, we'd need to move it, or `container.nf` to a directory on a lower level.

### Multiple Processes 

In any Bioinformatics pipeline, we'd like to combine these processses together. We can see this by running the following DSL-2 example: 

```
nextflow processes01/processes.nf
```

### Multiple inputs
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

### Asynchronous Processes 
* Note: that Nextflow processes define the execution of asynchronous tasks i.e. they are not executed one after another as they are written in the pipeline script as it would happen in a common imperative programming language.

In other words, the only thing keeping index from running before quantification is the fact that index output is used for quantification. 

We can see how this might work in `processes.nf` by commenting out: 

```
quantification(index.out.index, read_pairs_ch)
```

and un-commenting: 
```
// index_dir = "$baseDir/nextflow-training/data/index"
// quantification(index_dir, read_pairs_ch)
```
and then re-running the script with: 

```
nextflow processes01/processes.nf
```

you will notice both processes start simultaneously. This is because we simply pass an index file instead of waiting for one to be generate.

### Pulling Data from Other Sources
Thus far all of our data has come from this repositories local data directory. 

However, we can use the power of nextflow to pull from online sources as well. 

To do this, try runnning: 
```
nextflow processes02/processes_data.nf -c processes02/alternative_data.config
```
This script specifies another config that tells our script where some github hosted data is. We an then pull that data into our program and combine it with our local data. This is also another example of how config files stack as the `nextflow.config` in processes02 directory is also being used.

Notice how the output now identifies three runs of quantification. This features can be particularly powerful for pulling down test data or data into a variety of environments. 

### Running on LSF 

We can also use config files to specify profiles for different computing environments like the lsf. 

To do so, we first need to some setup.

First ssh into juno:
```
ssh -A -Y <user_name>@juno.mskcc.org
```

Re-clone the repo:
```
git clone https://github.com/mskcc-omics-workflows/nextflow-training.git
```

Load singularity and install nextflow 
```
module load singularity/3.7.1
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

We can then run: 
```
cd nextflow-training
nextflow processes03/processes_profiles.nf  -c processes03/alternative_profiles.config -profile lsf
```
Through using `alternative_profiles.config` alone, we were able to successfully submit our job to the cluster. This is because `alternative_profiles.config` specifies different computing profiles that set different configuration parameters. 

This same script can be run in your local environment with: 
```
nextflow processes03/processes_profiles.nf  -c processes03/alternative_profiles.config -profile docker
```


### Alignment, a CCI Example

Now it's time to try an example a little more specific to the needs of CCI. 

To start first install some additional test data by running: 
```
source install_data.sh
```

This will install the folder test_nucleo to our current directory. 

Next try running the following nextflow script:
```
nextflow run alignment/alignment.nf -entry test_alignment
```

Notice there is new parameter used named `-entry`, which specifies which workflow within the the nextflow script we'd like to enter. 

We can see these different workflows by taking a closer look at the `alignment.nf` script.

In addition to the different workflows there are a few other key elements to note. 

First are the include clauses: 
```
include { BWA_MEM      } from './modules/bwa_mem.nf'
include { PICARD_ADDORREPLACEREADGROUPS      } from './modules/picard_addorreplacereadgroups.nf'
```

These allow us to import nextflow processes for use in the workflow defined in `alignment.nf`. This proves useful as it is easy to see how a single nextflow script might get crowded without this functionality.

This script also shows off the power of Channel Operators. 

This proves to be powerful example of the uses of nextflow. And from here we have enought to write a bioinformatics pipeline. However, there are still a few weaknesses. 

For one, advanced.config is getting quite bloated. It is true we can specify multiple configs. However, it is not exactly clear what would be the best way to do this. In a similar note, it is also not clear how to best organize our modules and workflows. Or how testing should be handled in a formalized way? 

It seems we are missing some best practice standards and accompanying tools to help us navigate these obstacles. 

Of course, this is something anyone could be solved by an individual team applying good software enginneering practices. Alternatively, it turns out there is already an effort to do this via nf-core.

## Modularizing Alignment via Nf-Core

### Overview of Alignment as a Module

- Alignment in mskcc-omics-workflows/modules: https://github.com/mskcc-omics-workflows/modules/blob/feature/alignment/subworkflows/nf-core/alignment/main.nf 

### Adding a module
To get started, we can install the nf-core python package with the following: 

```
pip install nf-core==2.7.1
```

* Note please install 2.7.1 as this is what materials for this training have been made using. 2.9 is a recent release and I haven't vetted for compatibility.

Then, clone the mskcc-omics-workflows/modules and start a feature branch 
```
git clone https://github.com/mskcc-omics-workflows/modules.git
git checkout feature/<module_name>
```

Add a new module:
```
nf-core modules create <module_name>
```
Follow DSL2 module guidelines and Pull Request to develop. 

- See mergefastq example from Yu Hu for the basics: https://github.com/mskcc-omics-workflows/modules/pull/63

### Setting up a pipeline with mskcc-omics-workflows/modules

Creating a pipeline: 
```
git remote add origin https://github.com/mskcc-omics-workflows/pipeline_example.git
git branch -M main 
git push -u origin main
```

Installing a workflow: 

```
nf-core subworkflows --git-remote https://github.com/mskcc-omics-workflows/mskcc-modules.git --branch feature/alignment install alignment
```

## Additional Info 
Other topics of interest:
- Execution from online
- How to contribute to mskcc-omics-workflows/modules 


Pieces adapted from Sequera Labs Nextflow Training: https://github.com/seqeralabs 
