## HPC_T_Assembly
### Overview
 This script automates a high-performance computing (HPC) pipeline for RNA sequencing data analysis. It includes steps for data trimming, assembly, alignment, quantification, and clustering. The script is designed to handle paired-end RNA-seq data and generate several intermediate and final output files necessary for downstream analyses.

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

## Code
The main script is HPC_T_Assembly.py, it's a python code that takes the parameters set in the configuration site and generates the bash scripts for each step of the pipeline.

Main functions:
* getreqs(): This function retrieves the paths to required software tools like Trinity, CD-HIT, Salmon, etc., by searching for them in pre-defined locations.
 * mainhpc(threads): This is the main function of the script. It takes the number of available threads as input and performs the following:
   * Calls getreqs() to get software paths.
   * Checks if multiple species are present based on the HPC_T_Assembly_Data.txt file.
     * If multiple species are present, it generates separate folders for each species and creates an instance of the pipeline inside each folder to run them in parallel.
   * If a single species is present:
     * Processes the HPC_T_Assembly_Data.txt file to get read paths.
     * Generates the main assembly pipeline script pipeline.sh using these read paths and configuration files. This script includes steps for:
       * Fastq trimming with fastp.
       * Assembly with SPAdes.
       * Statistics with Trinity Stats.
       * Transcript clustering with CD-HIT.
       * Quantification with Salmon.
       * Transcript functional annotation with Corset.
       * Open reading frame (ORF) prediction with TransDecoder.
     * Generates scripts for individual analysis steps like Salmon indexing, Salmon quantification, Corset analysis, Hisat2 alignment, and TransDecoder prediction.
     * Generates a cleanup script cleanup.sh to organize the output files and remove temporary files.
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

