## HPC_T_Assembly
### Overview
 This script automates a high-performance computing (HPC) pipeline for RNA sequencing data analysis. It includes steps for data trimming, assembly, alignment, quantification, clustering, ORF predictions, and Transcript statistics. The script is designed to handle paired-end RNA-seq data and generate the intermediate and final output files necessary for downstream analyses. 

![IMG_20240910_001739_590](https://github.com/user-attachments/assets/dd8b4106-7200-4427-83ea-2c21c96859ac)


### Prerequisites

#### Web Interface
* Python3
* Flask

#### Main Code
* Slurm
* Python3

### USE
1. Install flask


    `pip install flask`


2. Launch web interface
   
    `python HPC_T_Assembly_Configuration.py`

3. Generate and Download Configuration files

4. Upload Config.zip, HPC_T_Assembly.py, and HPC_T_Assembly_Data.txt (path to left and right reads separated by ",")
5. Unzip Config.zip and transfer all configuration files into a folder called Config

    `unzip Config.zip ` 

    `mkdir Config` 

    `mv *.config.txt Config`


   
8. Launch HPC_T_Assembly.py

    `python HPC_T_Assembly.py`


## Multiple Species Workflow

![IMG_20240910_001739_623](https://github.com/user-attachments/assets/f81c55df-6169-4367-ab84-b55bc3d942d4)

### Execution for Multiple Species
* As shown in the previous graphic, when running the software for multiple species, a separate folder is created for each species, and the pipeline is executed in parallel within these folders.

### Execution for Single Species
The execution of the pipeline for a single specie proceeds as follows:

* Fastp is launched for each left, right reads pair. Generating the trimmed reads which will be used in the next step and an HTML file with the quality analysis before and after trimming which can be reviewed by the user.
* A single multithreaded process of RNASpades with as many threads as set in the HPC configuration is launched to assemble the trimmed reads.
* A multithreaded instance of the clusterer CD-HIT-EST is launched with the parameters configured in the web interface or with the default parameters -c 0.95 -n 10 -M 1600 taking as input the transcript generated by RNASpades and outputting another transcript file.
* Salmon builds the index and transcript quantifiers that will be used in Corset to generate its transcripts.
* Corset is launched as a single threaded process and its corresponding transcripts are generated.
* After all the transcripts are generated a single threaded process of trinity stats is launched per each transcript (Spades, CD-HIT, Corset).
* Hisat builds an index and then proceeds to align the reads. (Runs parallel to step 6)
* TransDecoder followed by TransDecoder.predict are launched to identify the candidate coding regions within the transcripted sequences.
* An initial cleanup script organizes the output into 4 folders: Transcripts, ORF, Statistics, Intermediate_Files.
* A final script eliminates the downloaded software that was required for the execution of the pipeline.

## Output Multiple species mode
* Specie1(Name)
  - Single Specie Output
* Specie2(Name)
  - Single Specie Output  
* .
* SpecieN(Name)
  - Single Specie Output

 ## Single Specie Output
 * Intermediate_Files
    - ASSEMBLY
    - Config
    - CorsetOutput
    - Data
    - FastP
    - Hisat2
    - Salmon
    - Scripts
    - slurmerr
    - slurmout 
 * ORF
 * Statistics
    - cdhitstats.txt
    - corsetstats.txt
    - spadestats.txt
 * Transcripts
    - transcripts_Corset.fasta
    - transcripts_cdhit.fasta
    - transcripts_rnaspades.fasta

## Code Breakdown
The main script is HPC_T_Assembly.py, it's a python code that takes the parameters set in the configuration site and generates the bash scripts for each step of the pipeline.

Main functions:
* getreqs(): This function retrieves the paths to required software tools like Trinity, CD-HIT, Salmon, etc., by searching for them in pre-defined locations.
 * mainhpc(threads): This is the main function of the script. It takes the number of available threads as input and performs the following:
   * Calls getreqs() to get software paths.
   * Checks if multiple species are present based on the HPC_T_Assembly_Data.txt file.
     * If multiple species are present, it generates separate folders for each species and creates an instance of the pipeline inside each folder to run them in parallel.
   * If a single species is present:
     * Processes the HPC_T_Assembly_Data.txt file to get read paths.
     * Generates the main assembly pipeline scripts using these read paths and configuration files. This script includes steps for:
       * Fastq trimming with fastp.
       * Assembly with SPAdes.
       * Statistics with Trinity Stats.
       * Transcript clustering with CD-HIT.
       * Quantification with Salmon.
       * Transcript functional annotation with Corset.
       * Open reading frame (ORF) prediction with TransDecoder.
 
 * cleanup(): This function creates a cleanup script that organizes output files, moves them to designated directories, and removes temporary files.
 * install_missing(name=None, instlist=None): This function is used to install required software tools from a predefined list based on their names or a provided list (in case some software didn't install correctly). It downloads and compiles the software from public repositories.
 * getsbatch(configf): This function retrieves the SLURM batch configuration (memory, time, etc.) from a specified configuration file.
 * getaccount(configf): This function retrieves the account to be used for running jobs on the HPC system from a configuration file.
 * getpartition(): This function retrieves the partition to be used for submitting jobs on the HPC system from a configuration file.
 * getcommand(configf): This function retrieves the actual command to be executed for a specific analysis step from a configuration file. This allows for customization of commands within the pipeline.
 * remove(): This function checks the sbatch.config.txt file to see if the remove_software step is to be executed.
Main block:
 * The code checks if any arguments are provided. (used internally in the code)
   * If no arguments are provided, it displays a menu asking the user to choose between immediate execution or manual configuration followed by execution.
   * If multispecie or Execute arguments are provided, it skips to the execution step.
 * Installs missing software tools if required (install_missing()).
 * Prepares the HPC_T_Assembly_Data.txt file if it doesn't exist by assuming paired-end reads based on filenames containing "_1" and "_2" located in a folder called Data.
 * Generates the main assembly pipeline script HPC_T_Assembly_Single.sh using a template and information from configuration files. This script includes steps to submit individual analysis steps using SLURM batch commands based on the defined configurations.
 * Generates a Processes.txt file summarizing the number of processes involved in each script.
 * Calls cleanup() to create the cleanup script.
 * Calls mainhpc(len(threadcounter(0))) to initiate the pipeline execution. 

