from os import getcwd as gc
from os.path import abspath, expandvars
from os import system
import yaml
from sys import argv
from time import time

reqd = {}
base = ['#SBATCH -N ', '#SBATCH -n ', '#SBATCH --mem=', '#SBATCH --account ', '#SBATCH --time ']
Multispecie = False
ExecuteNow = True if argv[1:2] in [["multispecie"], ["Execute"]] else False  # False if argv[1:2] == [] else True


def mainhpc(threads):
    reqd = getreqs()  # Get complete path to required software
    with open("HPC_T_Assembly_Data.txt") as f:
        allreads = f.read()
        l = allreads.split("#")
    if len(l) > 1:  # Determine if multiple species
        Multispecie = True
        left = ""
        right = ""
        for x in allreads.split("\n"):  # Split reads into left and right
            if "fastq" in x:
                left, right = left + x.split(",")[0], right + x.split(",")[1]
            else:
                left, right = left + x, right + x
            left, right = left + "\n", right + "\n"
        left, right = left.split("#"), right.split("#")

        species = []
        with open("HPC_T_Assembly_Multiple.sh", "w") as genscript:  # Generate script to run multiple species
            for specie, rreads in zip(left[1:], right[1:]):
                specie = specie.split("\n")
                specie_name = "_".join(specie[0].split())
                species.append(specie_name)
                lreads = [abspath(expandvars(x)) for x in specie[1:] if len(x) > 5]
                rreads = [abspath(expandvars(x)) for x in rreads.split("\n")[1:] if len(x) > 5]
                s(f"mkdir {specie_name}")
                s(f"cp -r Config {specie_name}")
                s(f"cp HPC_T_Assembly.py {specie_name}")
                with open(f"{specie_name}/HPC_T_Assembly_Data.txt", "w") as f:
                    f.write("\n".join(f"{x},{y}" for x, y in zip(lreads, rreads)))

                genscript.write(f"cd {specie_name}\npython HPC_T_Assembly.py multispecie\ncd ..\n")
        if ExecuteNow:  # If execute now --> run multispecie
            s("bash HPC_T_Assembly_Multiple.sh")
        exit()
    with open("HPC_T_Assembly_Data.txt") as f:  # For single specie
        rleft, rright = [], []
        for x in f.read().split("\n"):
            if len(x) > 5:
                rleft.append(x.split(",")[0])
                rright.append(x.split(",")[1])

    left, right = [], []
    for x, y in zip(rleft, rright):
        with open("HPC_T_Assembly_Data.txt", "w") as f:
            f.write(f"{abspath(expandvars(x))},{abspath(expandvars(y))}\n")
            left.append(abspath(expandvars(x)))
            right.append(abspath(expandvars(y)))

    threadsperrun = threads // len(left)
    with open("Config/fastp.config.txt", "r") as f:
        trim = f.read()

    t = [x.split("\n") for x in trim.split("#")]

    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]

    sbatchc = ""  # Generate sbatch config
    for x, y in zip(base, cvals):
        sbatchc += f'{x}{y}\n'

    for option in t[1][1:-1]:
        if len(option) > 1:
            sbatchc += f'#SBATCH {option}\n'

    command = " ".join(t[2][1:-1]) + " ".join(t[3][1:])
    left = [x for x in left if len(x) > 4]
    with open(f"pipeline.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(sbatchc)
        f.write(f'cd {gc()}\n')
        for i in range(len(left)):
            lefti, righti, leftio, rightio = left[i], right[i], left[i].split('.')[0], right[i].split('.')[0]
            f.write(command.format(**locals()) + "\n")

    # Assembly
    lreads = []
    rreads = []

    for i in range(len(left)):
        if len(left[i]) > 4:
            lreads.append(f"{left[i].split('.')[0]}_cleaned.fastq")
            rreads.append(f"{right[i].split('.')[0]}_cleaned.fastq")

    data = [
        {
            "orientation": "fr",
            "type": "paired-end",
            "right reads": rreads,
            "left reads": lreads
        }
    ]
    with open("fastq_files.yaml", "w") as f:
        yaml.safe_dump(data, f)

        # Spades
    # ../SPAdes-3.15.5-Linux
    spades = getcommand("Config/assembly.config.txt").format(**locals())
    # spades =f"{reqd['SPAdes']}/bin/rnaspades.py -t {threads} -o ASSEMBLY/ELA_spades_k_auto --dataset fastq_files.yaml --tmp-dir ASSEMBLY/TMP_SPADE/ --only-assembler"

    with open("assembly.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/assembly.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(spades)

    # Statistics
    # trinity = f"perl {reqd['trinityrnaseq']}/util/TrinityStats.pl ASSEMBLY/ELA_spades_k_auto/transcripts.fasta > spadestats.txt"
    trinity = getcommand("Config/trinitystats.config.txt").format(**locals())

    with open("statistics.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/trinitystats.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(trinity)
        f.write("\n" + " ".join(trinity.split()[0:2]) + " transcripts_cdhit.fasta > cdhitstats.txt")
        f.write("\n" + " ".join(trinity.split()[0:2]) + " transcripts_Corset.fasta > corsetstats.txt")

    # Clustering
    # CD-HIT-EST
    # cdhit = f"{reqd['cdhit']}/cd-hit-est -i transcripts.fasta -o transcripts.fasta -c 0.95 -n 10 -T {threads - 2}]"
    cdhit = getcommand("Config/cdhit.config.txt").format(**locals())

    with open("cdhit.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/cdhit.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(cdhit)

        # Corset
    # salmonpidx = f"{reqd['salmon']}/bin/salmon index --index Data/Salmon_Index --transcripts transcripts.fasta -p {threads}"
    salmonpidx = getcommand("Config/salmonidx.config.txt").format(**locals())
    salmon = []
    salmonpos = []
    for x, y in zip(lreads, rreads):
        ID = x.split("_cleaned.fastq")[0].split("/")[-1]
        # salmon.append(f'{reqd["salmon"]}/bin/salmon quant --index Data/Salmon_Index --libType A -1 {x} -2 {y} --dumpEq --output ELA_SALMON_{ID} &')
        salmon.append(getcommand("Config/salmon.config.txt").format(**locals()))
        # salmonpos.append(f'gunzip -k ELA_SALMON_{ID}/aux_info/eq_classes.txt.gz')
        salmonpos.append(getcommand("Config/salmonpos.config.txt").format(**locals()))

    with open("salmonidx.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/salmonidx.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(salmonpidx)

    with open("salmon.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/salmon.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(" &\n".join(salmon))

    with open("salmonpos.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/salmonpos.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(" &\n".join(salmonpos))

    with open("corset.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/Corset.config.txt"))
        f.write(f'cd {gc()}\n')
        # f.write(f'{reqd["corset"]}/corset corset -i salmon_eq_classes ELA_SALMON_{lreads[0].split("_cleaned.fastq")[0].split("/")[-1][:-1]}*/aux_info/eq_classes.txt -f true')
        ID = lreads[0].split("_cleaned.fastq")[0].split("/")[-1][:-1]
        f.write(getcommand("Config/Corset.config.txt").format(**locals()))

    with open("corset2transcript.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/Corset2transcript.config.txt"))
        # f.write(f"\npython {reqd['Corset-tools']}/fetchClusterSeqs.py -i transcripts.fasta -o transcripts_corset.fasta -c clusters.txt")
        f.write(getcommand("Config/Corset2transcript.config.txt").format(**locals()))

    # Final Statistics
    # Bowtie2
    # idxfasta = "transcripts_corset.fasta"
    # idxname = "Index"
    # index = f"{reqd['bowtie2']}/bowtie2-build {idxfasta} {idxname}"
    index = getcommand("Config/bowtie2index.config.txt").format(**locals())
    bowtie = []
    for x, y in zip(lreads, rreads):
        ID = x.split("_cleaned.fastq")[0].split("/")[-1].split("_")[0]

        bowtie.append(getcommand("Config/bowtie2.config.txt").format(**locals()))

    with open("bowtieindex.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/bowtie2index.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(index)
        # f.write(" &\n".join(bowtie))

    with open("bowtie2.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/bowtie2.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(" &\n".join(bowtie))

    # Busco
    busco = getcommand("Config/busco.config.txt").format(**locals())

    with open("busco.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/busco.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(busco)

        # ORF Predictions

    transdecoder_orf = getcommand("Config/transdecoder.config.txt").format(**locals())

    with open("transdecoder.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/transdecoder.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(transdecoder_orf)

    transdecoder_predict = getcommand("Config/transdecoder_predict.config.txt").format(**locals())

    with open("transdecoder_predict.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/transdecoder_predict.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(transdecoder_predict)

    # Remove Software
    with open("remove_software.sh", "w") as f:
        f.write("#!/bin/bash\n")

        if "Software" in ls():
            f.write("yes | rm -rf Software\n")
        else:
            f.write("yes | rm -rf ../Software\n")

    if argv[1:2] == ["Execute"] or ExecuteNow:
        s("bash HPC_T_Assembly_Single.sh")


def remove():
    with open("Config/sbatch.config.txt") as f:
        sbatchc = f.read()
    if "remove_software.sh" in sbatchc:
        return True
    else:
        return False


def getsbatch(configf):
    base = ['#SBATCH -N ', '#SBATCH -n ', '#SBATCH --mem=', '#SBATCH --account ', '#SBATCH --time ']
    with open(configf) as f:
        trim = f.read()
    t = [x.split("\n") for x in trim.split("#")]
    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]
    sbatchc = ""
    for x, y in zip(base, cvals):
        sbatchc += f'{x}{y}\n'
    for option in t[1][1:-1]:
        if len(option) > 1:
            sbatchc += f'#SBATCH {option}\n'
    return sbatchc


def getaccount(configf):
    with open(configf) as f:
        trim = f.read()
    t = [x.split("\n") for x in trim.split("#")]
    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]
    return cvals[3]


def getpartition():
    with open("Config/assembly.config.txt") as f:
        f = f.read()
    return [x for x in f.split("#")[1].split("\n") if x[0:2] == "-p"][0]


def getcommand(configf):
    with open(configf) as f:
        trim = f.read()
    t = [x.split("\n") for x in trim.split("#")]
    return " ".join(t[2][1:]) if len(t) < 4 else " ".join(t[2][1:-1]) + " ".join(t[3][1:]) if len(
        t[3]) > 1 else " ".join(t[2][1:-1])


def getreqs():
    ld = ls("/")
    reqs = """
    trinityrnaseq
    cdhit
    salmon
    corset-
    Corset-tools
    SPAdes
    fastp
    bowtie
    TransDecoder-
    busco
    """.split()
    mr = reqs
    reqs = {x: "" for x in reqs}
    try:
        ld = ls("Software")
        loc = "Software/"
    except:
        try:
            ld = ls("../Software")
            loc = "../Software/"
        except:
            install_missing()
            ld = ls("Software")
            loc = "Software/"
    for x in mr:
        for y in ld:
            if x in y and isdir(loc + y):
                reqs[x] = loc + y
                break
    return reqs


def cleanup():
    with open("Config/sbatch.config.txt") as f:
        sbatchc = f.read()
        if "Retries" not in sbatchc:
            retry = True
            lives=str(int(sbatchc.split("\n")[0].split(",")[-1]))
        else:
            retry = False

    with open("cleanup.sh", "w") as f:
        f.write("""#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH """ + getpartition() + """
#SBATCH --mem=12GB
#SBATCH --time 00:15:00
""" + f"#SBATCH --account {getaccount('Config/assembly.config.txt')}\n")
        if retry:
            f.write(f"remaining={lives}\n")
        else:
            f.write("remaining=0\n")
        f.write("""sl=0
fix() {
  if [[ $remaining -eq 0 ]]; then
    echo "Error: Could not fix the issue"
    exit 1
  fi
  line=${1:-0}
  line=$((line - sl)) 
  tail -n +"$((line + 1))" HPC_T_Assembly_Single.sh | sed '1s/--dependency=afterany:[^ ]*//' > HPC_T_Assembly_Single_tmp.sh
  mv HPC_T_Assembly_Single_tmp.sh HPC_T_Assembly_Single.sh
  sed -i -E "s/^sl=[0-9]+/sl=$line/" "cleanup.sh"
  sed -i -E "s/^remaining=[0-9]+/remaining=$((remaining - 1))/" "cleanup.sh"
  bash HPC_T_Assembly_Single.sh
}

# Read number of reads from Processes.txt
reads=$(sed -n '2p' Processes.txt | awk -F'|' '{print $2}')

# Verify fastp
cat fastp.err | grep CANCELLED > fastp.v
if grep -q "CANCELLED" fastp.v; then
  echo "FastP Verification Failed"
  fix 0
  exit 1
fi               
cat fastp.* | grep "Duplication rate:" > fastp.v
fastp_lines=$(wc -l < fastp.v)
if [ "$fastp_lines" -lt "$reads" ]; then
  echo "Fastp Verification Failed"
  fix 0
  exit 1
fi

# Verify assembly
cat assembly.* | grep "Assembling finished" > assembly.v
cat assembly.err | grep CANCELLED >> assembly.v
if grep -q "CANCELLED" assembly.v; then
  echo "Assembly Verification Failed"
  fix 1
  exit 1
fi  
if ! grep -q "Assembling finished." assembly.v; then
  echo "Assembly Verification Failed"
  fix 1
  exit 1
fi

# Verify CDHIT
cat cdhit.* | grep "writing new database
writing clustering information
program completed" > chdit.v
cat cdhit.err | grep CANCELLED >> cdhit.v
if grep -q "CANCELLED" cdhit.v; then
  echo "CDHIT Verification Failed"
  fix 2
  exit 1
fi  
if ! grep -q "writing new database
writing clustering information
program completed" cdhit.v; then
  echo "CDHIT Verification Failed"
  fix 2
  exit 1
fi

# Verify salmonidx
cat salmonidx.out | grep "Edges construction time:" > salmonidx.v
cat salmonidx.err | grep CANCELLED >> salmonidx.v
if grep -q "CANCELLED" salmonidx.v; then
  echo "SalmonIDX Verification Failed"
  fix 3
  exit 1
fi  
if ! grep -q "Edges construction time:" salmonidx.v; then
  echo "Salmon idx Verification Failed"
  fix 3
  exit 1
fi

# Verify salmon
cat salmon.err | grep CANCELLED > salmon.v
if grep -q "CANCELLED" salmon.v; then
  echo "Salmon Verification Failed"
  fix 4
  exit 1
fi  
cat salmon.* | grep "done writing equivalence class counts." > salmon.v
salmon_lines=$(wc -l < salmon.v)
if [ "$salmon_lines" -lt "$reads" ]; then
  echo "Salmon Verification Failed"
  fix 4
  exit 1
fi

# Verify corset
cat Corset.* | grep "Finished" > corset.v
cat Corset.err | grep CANCELLED >> corset.v
if grep -q "CANCELLED" corset.v; then
  echo "Corset Verification Failed"
  fix 6
  exit 1
fi  
if ! grep -q "Finished" corset.v; then
  echo "Corset Verification Failed"
  fix 6
  exit 1
fi
# Verify SalmonPos
cat salmonpos.err | grep "No such file or directory" > salmonpos.v
if grep -q "No such file or directory" salmonpos.v; then
    echo "SalmonPos Verification Failed"
    fix 5
    exit 1
fi
# Verify BowtieIDX
cat bowtieindex.err | grep CANCELLED > bowtie.v
if grep -q "CANCELLED" hisat.v; then
  echo "Hisat Verification Failed"
  fix 9
  exit 1
fi 
cat bowtie2index.err | grep "Renaming Index" > bowtie.v
if ! grep -q "Renaming Index" bowtie.v; then
    echo "Bowtie Index Verification Failed"
    fix 9
    exit 1
fi
  
# Verify Bowtie2
cat bowtie2.err | grep CANCELLED > bowtie2.v
if grep -q "CANCELLED" bowtie2.v; then
  echo "Bowtie2 Verification Failed"
  fix 11
  exit 1
fi  
cat bowtie2.* | grep "overall alignment rate" > bowtie.v
bowtie_lines=$(wc -l < bowtie.v)
if [ "$bowtie_lines" -lt "$reads" ]; then
  echo "Bowtie2 Verification Failed"
  fix 11
  exit 1
fi

#Verify Busco
cat busco.err | grep CANCELLED > busco.v
if grep -q "CANCELLED" busco.v; then
  echo "Busco Verification Failed"
  fix 10
  exit 1
fi
cat busco.* | grep "BUSCO analysis failed!" > busco.v
if grep -q "BUSCO analysis failed!" busco.v; then
    echo "Busco Verification Failed"
    fix 10
    exit 1
fi""")
        if remove():
            f.write("bash remove_software.sh\n")

        f.write("""mkdir Scripts
mv *.sh Scripts

mkdir slurmout
mv *.out slurmout

mkdir Salmon
mv -f RNA_SALMON* Salmon

mkdir bowtie2
mv *.sam bowtie2

rm r.txt
rm SRR*.txt

mkdir CorsetOutput
mv clusters.txt CorsetOutput
mv counts.txt CorsetOutput
cp transcripts_Corset.fasta CorsetOutput


mkdir ORF
mv -f transcripts_Corset.fasta.transdecoder* ORF


mkdir Data2
mv $(cat HPC_T_Assembly_Data.txt) Data2


mv -f Data/Salmon_Index Salmon

mv -f Data Fastp

mv -f Data2 Data

rm paths.txt

mkdir slurmerr
mv *.err slurmerr

mkdir Bowtie2Output 
mv log_*.txt Bowtie2Output
mkdir Transcripts
mv transcripts_cdhit.fasta Transcripts
cp transcripts_Corset.fasta Transcripts
mv transcripts.fasta Transcripts/transcripts_rnaspades.fasta  
mkdir Intermediate_Files
mv * Intermediate_Files
mv Intermediate_Files/Transcripts .
mv Intermediate_Files/ORF .   
mkdir Statistics
mv Intermediate_Files/*stats.txt Statistics 
""")



def install_missing(name=None, instlist=None):  # Install required software

    installation_instructions = {
        'trinityrnaseq': '\ngit clone https://github.com/trinityrnaseq/trinityrnaseq.git\ncd trinityrnaseq\nmake\nmake plugins\ncd ..\n',
        'bowtie': 'wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-sra-linux-x86_64.zip -O bowtiebin.zip\nunzip bowtiebin.zip\nmv bowtie2-2.5.4-sra-linux-x86_64 bowtie2-2.5.4\n',
        # 'bowtie': '\nwget https://sourceforge.net/projects/bowtie-bio/files/latest/download -O bowtie.zip\nunzip bowtie.zip\nmkdir bowtie\nmv bowtie*/* bowtie\nrm -rf bowtie2-*\ncd bowtie\ncmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .\n',
        'salmon': '\nwget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz\ntar -zxvf salmon-1.10.0_linux_x86_64.tar.gz\n',
        'samtools': '\ngit clone https://github.com/samtools/samtools.git\ncd samtools\n./configure\nmake\nmake install\ncd ..\n',
        'cdhit': '\ngit clone https://github.com/weizhongli/cdhit.git\ncd cdhit\nmake\ncd cd-hit-auxtools\nmake\ncd ..\ncd ..\n',
        'corset-': '\nwget https://github.com/Oshlack/Corset/releases/download/version-1.09/corset-1.09-linux64.tar.gz\ntar -zxvf corset-1.09-linux64.tar.gz\n',
        'Corset-tools': '\ngit clone https://github.com/Adamtaranto/Corset-tools.git\n',
        'SPAdes': '\nwget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz\ntar -zxvf SPAdes-4.0.0-Linux.tar.gz\n',
        'fastp': '\nwget http://opengene.org/fastp/fastp\nchmod +x fastp\nmkdir fastp1\nmv fastp fastp1\nmv fastp1 fastp\n',
        'busco': '\ngit clone https://gitlab.com/ezlab/busco.git\ncd busco/\npython -m pip install .\npip install pandas\npip install requests\npip install biopython\nwget https://sourceforge.net/projects/bbmap/files/latest/download\ntar -zxvf download\nexport PATH=$(pwd)/bbmap:$PATH\nwget https://mmseqs.com/metaeuk/metaeuk-linux-avx2.tar.gz; tar xzvf metaeuk-linux-avx2.tar.gz; export PATH=$(pwd)/metaeuk/bin/:$PATH\n',
        'TransDecoder-': '\ncurl -L https://cpanmin.us | perl - App::cpanminus\ncpanm install DB_Filea\ncpanm install URI::Escape\nwget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.zip\nunzip TransDecoder-v5.7.1.zip\n'}
    mr = list(installation_instructions.keys()) if instlist == None else instlist
    if name:
        s("bash Software/" + name)
        return
    if argv[1:2] == [] or name == False:  # If used to install all the software
        for tool in mr:
            with open(f"{tool}.sh", "w") as f:
                f.write(installation_instructions[tool])
                f.write(f"mv {tool}.sh Scripts\n")
        with open("install.sh", "w") as f:
            f.write("mkdir Software\n")
            f.write("cd Software\n")
            f.write("mkdir Scripts\n")
            f.write("mv ../*.sh .\n")
            mr = ["./" + x for x in mr]
            f.write(".sh &\n".join(mr))
            f.write(".sh\n")

        system("chmod +x *.sh")
        system("./install.sh")

        ct = time()  # Wait for installation to complete
        while not all(getreqs().values()):
            if (time() - ct) % 60 == 0:
                for missingreq in getreqs().items():
                    if missingreq[1] == "":
                        install_missing(missingreq[0] + ".sh")
            if (time() - ct) % 300 == 0:
                install_missing(False, [x[0] for x in getreqs().items() if x[1] == ""])

            print("Waiting for installation to complete")
            sleep(10)

        if ExecuteNow:
            s("python HPC_T_Assembly.py Execute")  # Execute after installation
        else:
            s("python HPC_T_Assembly.py False")  # Generate scripts and exit
        exit()


from os import system as s
from os import listdir as ls
from os.path import isdir
from os import sched_getaffinity as threadcounter
import yaml
from datetime import datetime
from time import sleep

if __name__ == "__main__":
    if argv[1:2] == []:
        print("=== Main Menu ===")
        print("1. Execute now\n2. Edit Configuration Batch and Execute Manually")
        if input(": ") == "1":
            ExecuteNow = True
    if argv[1:2] == []:
        install_missing()
    if "HPC_T_Assembly_Data.txt" not in ls():
        s("ls Data/*.fastq > r.txt")
        with open("r.txt") as f:
            ld = f.read()

        ld = ld.split()

        left = []
        right = []
        for x in ld:
            if "_1" in x:
                left.append(x)
            elif "_2" in x:
                right.append(x)

        with open("HPC_T_Assembly_Data.txt", "w") as f:
            f.write("\n".join(f"{x},{y}" for x, y in zip(left, right)))

    with open("Config/sbatch.config.txt") as f:
        sbatchc = f.read()

    sbatchcmd = ""
    for line in sbatchc.split("\n")[1:]:
        if len(line) > 5:
            cline = line.split()
            if len(cline) == 2:
                sbatchcmd += f'{cline[0].split(".")[0]}=$(sbatch --parsable --mem={cline[1]} {cline[0]})\n'
            else:
                sbatchcmd += f'{cline[0].split(".")[0]}=$(sbatch --parsable --dependency=afterany:${":".join(['{' + x.split(".")[0] + '}' for x in cline[1:-1]])} --mem={cline[-1]} {cline[0]})\n'

    with open("HPC_T_Assembly_Single.sh", "w") as f:
        f.write(sbatchcmd)

    with open("HPC_T_Assembly_Data.txt") as f:
        left = f.read().split("\n")

    lenleft = len(left)

    with open("Processes.txt", "w") as f:
        f.write("Script | Number of Processes\n")
        f.write(
            f"pipeline.sh | {lenleft}\nassembly.sh | 1\nsalmonidx.sh | 1\nsalmon.sh | {lenleft}\nsalmonpos.sh | {lenleft}\ncorset.sh | 1\ncorset2transcript.sh | 1\nbowtieindex.sh | 1\nbowtie2.sh | {lenleft}\ntrasdecoder.sh | 1\ntransdecoder_predict.sh | 1\n")
    cleanup()
    from os import sched_getaffinity as threadcounter

    mainhpc(len(threadcounter(0)))
