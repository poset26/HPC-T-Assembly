from flask import Flask, request, render_template, send_file, make_response
import os

#app = Flask(__name__)
app = Flask(__name__, static_folder='static')

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/main')
def main():
    return render_template('main.html')

@app.route('/generate', methods=['POST'])
def generate():  
    sbatchconfig = {"nodes":request.form["nodes"],"threads":request.form["threads"], "memory":request.form["memory"], "account":request.form["account"], "time":request.form["time"], "partition":request.form["p_flag"]} 
    sbatch = ['pipeline.sh  ', 'assembly.sh pipeline.sh ', 'cdhit.sh assembly.sh ', 'salmonidx.sh cdhit.sh ', 'salmon.sh salmonidx.sh ', 'salmonpos.sh salmon.sh ', 'corset.sh salmonpos.sh ', 'corset2transcript.sh corset.sh ', 'statistics.sh corset2transcript.sh ', 'hisatindex.sh corset2transcript.sh ', 'hisat2.sh hisatindex.sh ', 'cleanup.sh hisat2.sh 12g']
    sbatcht = sbatch[:-1] + ['transdecoder.sh hisat2.sh ', 'transdecoder_predict.sh transdecoder.sh ', 'cleanup.sh transdecoder_predict.sh 12g']
    stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "hisat2index", "hisat2"] #Use amount of threads specified by the user
    singlecore = ["Corset", "Corset2transcript","trinitystats","transdecoder", "transdecoder_predict"] #Use only 1 thread
    addconfigs = ["fastp","assembly","cdhit","Corset","hisat2","salmon"]  
    sbatch = sbatcht if request.form.get("ORFs") == "true" else sbatch #Add ORF prediction to the pipeline if requested
    for file_name in stdfiles:
        with open(f"{file_name}.config.txt", 'w') as file:
            file.write(f"Nodes: {sbatchconfig['nodes']}\n")
            file.write(f"Threads: {sbatchconfig['threads']}\n")
            file.write(f"Memory: {sbatchconfig['memory']}\n")
            file.write(f"Account: {sbatchconfig['account']}\n")
            file.write(f"Time: {sbatchconfig['time']}\n")
            file.write("#Other Sbatch configs\n")
            file.write(f"-p {sbatchconfig['partition']}\n")
            file.write(f"-o {file_name}.out\n-e {file_name}.err\n")
            file.write(f"# {file_name}\n")
            file.write(defaultconfigs[file_name])
            if file_name in addconfigs:
                file.write(f"\n{request.form[f'{file_name}_additional_configurations']}")
            
    for file_name in singlecore:
        with open(f"{file_name}.config.txt", 'w') as file:
            file.write(f"Nodes: 1\n")
            file.write(f"Threads: 1\n")
            file.write(f"Memory: {sbatchconfig['memory']}\n")
            file.write(f"Account: {sbatchconfig['account']}\n")
            file.write(f"Time: {sbatchconfig['time']}\n")
            file.write("#Other Sbatch configs\n")
            file.write(f"-p {sbatchconfig['partition']}\n")
            file.write(f"-o {file_name}.out\n-e {file_name}.err\n")
            file.write(f"# {file_name}\n")
            file.write(defaultconfigs[file_name])
            if file_name in addconfigs:
                file.write(f"\n{request.form[f'{file_name}_additional_configurations']}")
    with open("sbatch.config.txt", "w") as f:
        f.write("# Script, Dependencies, Memory \n")
        f.write(f"{sbatchconfig['memory'][:-2]}g\n".join(sbatch))
        if request.form.get("Delete") == "true":
            f.write("\nremove_software.sh cleanup.sh 12g\n")
    filelist = [x + ".config.txt" for x in stdfiles + singlecore] + ["sbatch.config.txt", "HPC_T_Assembly.py"]
    with open("HPC_T_Assembly.py","wb") as f:
        f.write(hpcpipeline)
    import zipfile

    with zipfile.ZipFile('Config.zip', 'w') as zipMe:
        for file in filelist:
            zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)
    
    return send_file(f"{os.getcwd()}/Config.zip", as_attachment=True)

hpcpipeline = b'from os import getcwd as gc\nfrom os.path import abspath, expandvars\nfrom os import system\nimport yaml\nfrom sys import argv\nfrom time import time\nreqd = {}\nbase = [\'#SBATCH -N \', \'#SBATCH -n \', \'#SBATCH --mem=\', \'#SBATCH --account \', \'#SBATCH --time \']\nMultispecie = False\nExecuteNow = True if argv[1:2] in [["multispecie"],["Execute"]] else False #False if argv[1:2] == [] else True\n\ndef mainhpc(threads):\n    reqd = getreqs()  # Get complete path to required software\n    with open("HPC_T_Assembly_Data.txt") as f:\n        allreads = f.read()\n        l = allreads.split("#")\n    if len(l) > 1:  # Determine if multiple species\n        Multispecie = True\n        left = ""\n        right = ""\n        for x in allreads.split("\\n"): # Split reads into left and right\n            if "fastq" in x:\n                left, right = left + x.split(",")[0], right + x.split(",")[1]\n            else:\n                left, right = left + x, right + x\n            left, right = left + "\\n", right + "\\n"\n        left, right = left.split("#"), right.split("#")\n\n\n        species = []\n        with open("HPC_T_Assembly_Multiple.sh", "w") as genscript: # Generate script to run multiple species\n            for specie, rreads in zip(left[1:], right[1:]):\n                specie = specie.split("\\n")\n                specie_name = "_".join(specie[0].split())\n                species.append(specie_name)\n                lreads = [abspath(expandvars(x)) for x in specie[1:] if len(x) > 5]\n                rreads = [abspath(expandvars(x)) for x in rreads.split("\\n")[1:] if len(x) > 5]\n                s(f"mkdir {specie_name}")\n                s(f"cp -r Config {specie_name}")\n                s(f"cp HPC_T_Assembly.py {specie_name}")\n                with open(f"{specie_name}/HPC_T_Assembly_Data.txt", "w") as f:\n                    f.write("\\n".join(f"{x},{y}" for x,y in zip(lreads,rreads)))\n\n                genscript.write(f"cd {specie_name}\\npython HPC_T_Assembly.py multispecie\\ncd ..\\n")\n        if ExecuteNow: # If execute now --> run multispecie\n            s("bash HPC_T_Assembly_Multiple.sh")\n        exit()\n    with open("HPC_T_Assembly_Data.txt") as f: # For single specie\n        rleft, rright = [],[]\n        for x in f.read().split("\\n"):\n            if len(x) > 5:\n                rleft.append(x.split(",")[0])\n                rright.append(x.split(",")[1])\n    \n    left, right = [],[]\n    for x,y in zip(rleft,rright):\n        with open("HPC_T_Assembly_Data.txt", "w") as f:\n            f.write(f"{abspath(expandvars(x))},{abspath(expandvars(y))}\\n")\n            left.append(abspath(expandvars(x)))\n            right.append(abspath(expandvars(y)))\n \n\n    threadsperrun = threads // len(left)\n    with open("Config/fastp.config.txt", "r") as f:\n        trim = f.read()\n\n    t = [x.split("\\n") for x in trim.split("#")]\n\n    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]\n\n    sbatchc = "" # Generate sbatch config\n    for x, y in zip(base, cvals):\n        sbatchc += f\'{x}{y}\\n\'\n\n    for option in t[1][1:-1]:\n        if len(option) > 1:\n            sbatchc += f\'#SBATCH {option}\\n\'\n\n    command = " ".join(t[2][1:-1]) + " ".join(t[3][1:])\n    left = [x for x in left if len(x) > 4]\n    with open(f"pipeline.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(sbatchc)\n        f.write(f\'cd {gc()}\\n\')\n        for i in range(len(left)):\n            lefti, righti, leftio, rightio = left[i], right[i], left[i].split(\'.\')[0], right[i].split(\'.\')[0]\n            f.write(command.format(**locals()) + "\\n")\n\n    # Assembly\n    lreads = []\n    rreads = []\n\n    for i in range(len(left)):\n        if len(left[i]) > 4:\n            lreads.append(f"{left[i].split(\'.\')[0]}_cleaned.fastq")\n            rreads.append(f"{right[i].split(\'.\')[0]}_cleaned.fastq")\n\n    data = [\n        {\n            "orientation": "fr",\n            "type": "paired-end",\n            "right reads": rreads,\n            "left reads": lreads\n        }\n    ]\n    with open("fastq_files.yaml", "w") as f:\n        yaml.safe_dump(data, f)\n\n        # Spades\n    # ../SPAdes-3.15.5-Linux\n    spades = getcommand("Config/assembly.config.txt").format(**locals())\n    # spades =f"{reqd[\'SPAdes\']}/bin/rnaspades.py -t {threads} -o ASSEMBLY/ELA_spades_k_auto --dataset fastq_files.yaml --tmp-dir ASSEMBLY/TMP_SPADE/ --only-assembler"\n\n    with open("assembly.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/assembly.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(spades)\n\n    # Statistics\n    # trinity = f"perl {reqd[\'trinityrnaseq\']}/util/TrinityStats.pl ASSEMBLY/ELA_spades_k_auto/transcripts.fasta > spadestats.txt"\n    trinity = getcommand("Config/trinitystats.config.txt").format(**locals())\n\n    with open("statistics.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/trinitystats.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(trinity)\n        f.write("\\n" + " ".join(trinity.split()[0:2]) + " transcripts_cdhit.fasta > cdhitstats.txt")\n        f.write("\\n" + " ".join(trinity.split()[0:2]) + " transcripts_Corset.fasta > corsetstats.txt")\n    \n    \n    \n    \n    # Clustering\n    # CD-HIT-EST\n    # cdhit = f"{reqd[\'cdhit\']}/cd-hit-est -i transcripts.fasta -o transcripts.fasta -c 0.95 -n 10 -T {threads - 2}]"\n    cdhit = getcommand("Config/cdhit.config.txt").format(**locals())\n\n    with open("cdhit.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/cdhit.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(cdhit)\n\n        # Corset\n    # salmonpidx = f"{reqd[\'salmon\']}/bin/salmon index --index Data/Salmon_Index --transcripts transcripts.fasta -p {threads}"\n    salmonpidx = getcommand("Config/salmonidx.config.txt").format(**locals())\n    salmon = []\n    salmonpos = []\n    for x, y in zip(lreads, rreads):\n        ID = x.split("_cleaned.fastq")[0].split("/")[-1]\n        # salmon.append(f\'{reqd["salmon"]}/bin/salmon quant --index Data/Salmon_Index --libType A -1 {x} -2 {y} --dumpEq --output ELA_SALMON_{ID} &\')\n        salmon.append(getcommand("Config/salmon.config.txt").format(**locals()))\n        # salmonpos.append(f\'gunzip -k ELA_SALMON_{ID}/aux_info/eq_classes.txt.gz\')\n        salmonpos.append(getcommand("Config/salmonpos.config.txt").format(**locals()))\n\n    with open("salmonidx.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/salmonidx.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(salmonpidx)\n\n    with open("salmon.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/salmon.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(" &\\n".join(salmon))\n\n\n    with open("salmonpos.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/salmonpos.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(" &\\n".join(salmonpos))\n\n    with open("corset.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/Corset.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        # f.write(f\'{reqd["corset"]}/corset corset -i salmon_eq_classes ELA_SALMON_{lreads[0].split("_cleaned.fastq")[0].split("/")[-1][:-1]}*/aux_info/eq_classes.txt -f true\')\n        ID = lreads[0].split("_cleaned.fastq")[0].split("/")[-1][:-1]\n        f.write(getcommand("Config/Corset.config.txt").format(**locals()))\n\n    with open("corset2transcript.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/Corset2transcript.config.txt"))\n        # f.write(f"\\npython {reqd[\'Corset-tools\']}/fetchClusterSeqs.py -i transcripts.fasta -o transcripts_corset.fasta -c clusters.txt")\n        f.write(getcommand("Config/Corset2transcript.config.txt").format(**locals()))\n\n    # Final Statistics\n    # HiSat2\n    # idxfasta = "transcripts_corset.fasta"\n    # idxname = "Index"\n    # index = f"{reqd[\'hisat2\']}/hisat2-build {idxfasta} {idxname}"\n    index = getcommand("Config/hisat2index.config.txt").format(**locals())\n    hisat = []\n    for x, y in zip(lreads, rreads):\n        ID = x.split("_cleaned.fastq")[0].split("/")[-1].split("_")[0]\n        # hisat.append(f\'{reqd["hisat2"]}/hisat2 -p {threadsperrun} --dta -q -x Index -1 {x} -2 {y} -S {x.split("_cleaned.fastq")[0].split("/")[-1].split("_")[0]}.sam > {x.split("_cleaned.fastq")[0].split("/")[-1].split("_")[0]}.txt\')\n        hisat.append(getcommand("Config/hisat2.config.txt").format(**locals()))\n\n    with open("hisatindex.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/hisat2index.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(index)\n        # f.write(" &\\n".join(hisat))\n\n    with open("hisat2.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/hisat2.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(" &\\n".join(hisat))\n\n    # ORF Predictions\n\n    transdecoder_orf = getcommand("Config/transdecoder.config.txt").format(**locals())\n\n    with open("transdecoder.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/transdecoder.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(transdecoder_orf)\n\n    transdecoder_predict = getcommand("Config/transdecoder_predict.config.txt").format(**locals())\n\n    with open("transdecoder_predict.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        f.write(getsbatch("Config/transdecoder_predict.config.txt"))\n        f.write(f\'cd {gc()}\\n\')\n        f.write(transdecoder_predict)\n        \n    # Remove Software\n    with open("remove_software.sh", "w") as f:\n        f.write("#!/bin/bash\\n")\n        \n        if "Software" in ls():\n            f.write("yes | rm -rf Software\\n")\n        else:\n            f.write("yes | rm -rf ../Software\\n")\n    \n    if argv[1:2] == ["Execute"] or ExecuteNow:\n        s("bash HPC_T_Assembly_Single.sh")\n    \n\ndef remove():\n    with open("Config/sbatch.config.txt") as f:\n        sbatchc = f.read()\n    if "remove_software.sh" in sbatchc:\n        return True\n    else:\n        return False\ndef getsbatch(configf):\n    base = [\'#SBATCH -N \', \'#SBATCH -n \', \'#SBATCH --mem=\', \'#SBATCH --account \', \'#SBATCH --time \']\n    with open(configf) as f:\n        trim = f.read()\n    t = [x.split("\\n") for x in trim.split("#")]\n    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]\n    sbatchc = ""\n    for x, y in zip(base, cvals):\n        sbatchc += f\'{x}{y}\\n\'\n    for option in t[1][1:-1]:\n        if len(option) > 1:\n            sbatchc += f\'#SBATCH {option}\\n\'\n    return sbatchc\n\ndef getaccount(configf):\n    with open(configf) as f:\n        trim = f.read()\n    t = [x.split("\\n") for x in trim.split("#")]\n    cvals = [x.split(":")[1].strip(" ") for x in t[0][:-2]] + [t[0][-2].split("Time: ")[1]]\n    return cvals[3]\n\ndef getpartition():\n    with open("Config/assembly.config.txt") as f:\n        f = f.read()\n    return [x for x in f.split("#")[1].split("\\n") if x[0:2] == "-p"][0]\n\ndef getcommand(configf):\n    with open(configf) as f:\n        trim = f.read()\n    t = [x.split("\\n") for x in trim.split("#")]\n    return " ".join(t[2][1:]) if len(t) < 4 else " ".join(t[2][1:-1]) + " ".join(t[3][1:]) if len(t[3]) > 1 else " ".join(t[2][1:-1])\n\n\ndef getreqs():\n    ld = ls("/")\n    reqs = """\n    trinityrnaseq\n    cdhit\n    salmon\n    corset-\n    Corset-tools\n    SPAdes\n    fastp\n    hisat2\n    TransDecoder-\n    """.split()\n    mr = reqs\n    reqs = {x: "" for x in reqs}\n    try:\n        ld = ls("Software")\n        loc = "Software/"\n    except:\n        try:\n            ld = ls("../Software")\n            loc = "../Software/"\n        except:\n            install_missing()\n            ld = ls("Software")\n            loc = "Software/"\n    for x in mr:\n        for y in ld:\n            if x in y and isdir(loc + y):\n                reqs[x] = loc + y\n                break\n    return reqs\n\n\ndef cleanup():\n    with open("cleanup.sh", "w") as f:\n        f.write("""#!/bin/bash\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH """ + getpartition() + """\n#SBATCH --mem=12GB\n#SBATCH --time 00:15:00\n""" + f"#SBATCH --account {getaccount(\'Config/assembly.config.txt\')}\\n")\n        if remove():\n            f.write("bash remove_software.sh\\n")\n        f.write("""\nmkdir Scripts\nmv *.sh Scripts\n\nmkdir slurmout\nmv *.out slurmout\n\nmkdir Salmon\nmv -f RNA_SALMON* Salmon\n\nmkdir Hisat2\nmv *.sam Hisat2\n\nrm r.txt\nrm SRR*.txt\n\nmkdir CorsetOutput\nmv clusters.txt CorsetOutput\nmv counts.txt CorsetOutput\ncp transcripts_Corset.fasta CorsetOutput\n\n\nmkdir ORF\nmv -f transcripts_Corset.fasta.transdecoder* ORF\n\n\nmkdir Data2\nmv $(cat HPC_T_Assembly_Data.txt) Data2\n\n\nmv -f Data/Salmon_Index Salmon\n\nmv -f Data Fastp\n\nmv -f Data2 Data\n\nrm paths.txt\n\nmkdir slurmerr\nmv *.err slurmerr\n\nmkdir Transcripts\nmv transcripts_cdhit.fasta Transcripts\ncp transcripts_Corset.fasta Transcripts\nmv transcripts.fasta Transcripts/transcripts_rnaspades.fasta  \nmkdir Intermediate_Files\nmv * Intermediate_Files\nmv Intermediate_Files/Transcripts .\nmv Intermediate_Files/ORF .   \nmkdir Statistics\nmv Intermediate_Files/*stats.txt Statistics \n""")\n\n\ndef install_missing(name = None, instlist = None): # Install required software\n\n    installation_instructions = {\n        \'trinityrnaseq\': \'\\ngit clone https://github.com/trinityrnaseq/trinityrnaseq.git\\ncd trinityrnaseq\\nmake\\nmake plugins\\ncd ..\\n\',\n        \'salmon\': \'\\nwget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz\\ntar -zxvf salmon-1.10.0_linux_x86_64.tar.gz\\n\',\n        \'samtools\': \'\\ngit clone https://github.com/samtools/samtools.git\\ncd samtools\\n./configure\\nmake\\nmake install\\ncd ..\\n\',\n        \'cdhit\': \'\\ngit clone https://github.com/weizhongli/cdhit.git\\ncd cdhit\\nmake\\ncd cd-hit-auxtools\\nmake\\ncd ..\\ncd ..\\n\',\n        \'corset-\': \'\\nwget https://github.com/Oshlack/Corset/releases/download/version-1.09/corset-1.09-linux64.tar.gz\\ntar -zxvf corset-1.09-linux64.tar.gz\\n\',\n        \'Corset-tools\': \'\\ngit clone https://github.com/Adamtaranto/Corset-tools.git\\n\',\n        \'SPAdes\': \'\\nwget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz\\ntar -zxvf SPAdes-4.0.0-Linux.tar.gz\\n\',\n        \'fastp\': \'\\nwget http://opengene.org/fastp/fastp\\nchmod +x fastp\\nmkdir fastp1\\nmv fastp fastp1\\nmv fastp1 fastp\\n\',\n        \'hisat2\': \'\\nwget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download\\nunzip download\\n\\n\',\n        \'TransDecoder-\': \'\\ncurl -L https://cpanmin.us | perl - App::cpanminus\\ncpanm install DB_Filea\\ncpanm install URI::Escape\\nwget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.zip\\nunzip TransDecoder-v5.7.1.zip\\n\'}\n    mr = list(installation_instructions.keys()) if instlist == None else instlist\n    if name:\n        s("bash Software/" + name)\n        return\n    if argv[1:2] == [] or name == False: # If used to install all the software\n        for tool in mr:\n            with open(f"{tool}.sh", "w") as f:\n                f.write(installation_instructions[tool])\n                f.write(f"mv {tool}.sh Scripts\\n")\n        with open("install.sh", "w") as f:\n            f.write("mkdir Software\\n")\n            f.write("cd Software\\n")\n            f.write("mkdir Scripts\\n")\n            f.write("mv ../*.sh .\\n")\n            mr = ["./" + x for x in mr]\n            f.write(".sh &\\n".join(mr))\n            f.write(".sh\\n")\n            \n        system("chmod +x *.sh")\n        system("./install.sh")\n\n        ct = time() # Wait for installation to complete\n        while not all(getreqs().values()):\n            if (time() - ct)%60 == 0:\n                for missingreq in getreqs().items():\n                    if x[1] == "":\n                        install_missing(x[0] + ".sh")\n            if (time() - ct)%300 == 0:\n                install_missing(False,[x[0] for x in getreqs().items() if x[1] == ""])\n                \n            print("Waiting for installation to complete")\n            sleep(10)\n            \n        \n        if ExecuteNow:\n            s("python HPC_T_Assembly.py Execute") # Execute after installation\n        else:\n            s("python HPC_T_Assembly.py False") # Generate scripts and exit\n        exit()\n\nfrom os import system as s\nfrom os import listdir as ls\nfrom os.path import isdir\nfrom os import sched_getaffinity as threadcounter\nimport yaml\nfrom datetime import datetime\nfrom time import sleep\n\n\nif __name__ == "__main__":\n    if argv[1:2] == []:\n        print("=== Main Menu ===")\n        print("1. Execute now\\n2. Edit Configuration Batch and Execute Manually")\n        if input(": ") == "1":\n            ExecuteNow = True\n    if argv[1:2] == []:\n        install_missing()\n    if "HPC_T_Assembly_Data.txt" not in ls():\n        s("ls Data/*.fastq > r.txt")\n        with open("r.txt") as f:\n            ld = f.read()\n\n        ld = ld.split()\n\n        left = []\n        right = []\n        for x in ld:\n            if "_1" in x:\n                left.append(x)\n            elif "_2" in x:\n                right.append(x)\n\n        with open("HPC_T_Assembly_Data.txt", "w") as f:\n            f.write("\\n".join(f"{x},{y}" for x,y in zip(left,right)))\n\n    with open("Config/sbatch.config.txt") as f:\n        sbatchc = f.read()\n\n    sbatchcmd = ""\n    for line in sbatchc.split("\\n")[1:]:\n        if len(line) > 5:\n            cline = line.split()\n            if len(cline) == 2:\n                sbatchcmd += f\'{cline[0].split(".")[0]}=$(sbatch --parsable --mem={cline[1]} {cline[0]})\\n\'\n            else:\n                sbatchcmd += f\'{cline[0].split(".")[0]}=$(sbatch --parsable --dependency=afterany:${":".join([\'{\' + x.split(".")[0] + \'}\' for x in cline[1:-1]])} --mem={cline[-1]} {cline[0]})\\n\'\n\n    with open("HPC_T_Assembly_Single.sh", "w") as f:\n        f.write(sbatchcmd)\n\n    with open("HPC_T_Assembly_Data.txt") as f:\n        left = f.read().split("\\n")\n\n    lenleft = len(left)\n\n    with open("Processes.txt", "w") as f:\n        f.write("Script | Number of Processes\\n")\n        f.write(\n            f"pipeline.sh | {lenleft}\\nassembly.sh | 1\\nsalmonidx.sh | 1\\nsalmon.sh | {lenleft}\\nsalmonpos.sh | {lenleft}\\ncorset.sh | 1\\ncorset2transcript.sh | 1\\nhisatindex.sh | 1\\nhisat2.sh | {lenleft}\\ntrasdecoder.sh | 1\\ntransdecoder_predict.sh | 1\\n")\n    cleanup()\n    from os import sched_getaffinity as threadcounter\n\n    mainhpc(len(threadcounter(0)))\n'

# file configurations
stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "hisat2index", "hisat2"]
singlecore = ["Corset", "Corset2transcript","trinitystats","transdecoder", "transdecoder_predict"]
defaultconfigs = {
singlecore[0]:"""{reqd[corset-]}/corset corset
-i salmon_eq_classes RNA_SALMON_*/aux_info/eq_classes.txt -f true""",
singlecore[1]:"""python {reqd[Corset-tools]}/fetchClusterSeqs.py 
-i transcripts_cdhit.fasta
-o transcripts_Corset.fasta
-c clusters.txt""",
singlecore[2]:"""perl {reqd[trinityrnaseq]}/util/TrinityStats.pl
transcripts.fasta
> spadestats.txt""",
singlecore[3]:"""{reqd[TransDecoder-]}/TransDecoder.LongOrfs 
-t transcripts_Corset.fasta""",
singlecore[4]:"""{reqd[TransDecoder-]}/TransDecoder.Predict
-t transcripts_Corset.fasta""",
stdfiles[0]:"""{reqd[fastp]}/fastp
-i {lefti}
-I {righti}
-o {leftio}_cleaned.fastq 
-O {rightio}_cleaned.fastq 
-w 48
-j {leftio}_fastp.json
-h {leftio}_fastp.html""",
stdfiles[1]:"""{reqd[SPAdes]}/bin/rnaspades.py
-t 48
-o ASSEMBLY/RNA-Spades""",
stdfiles[2]:"""{reqd[cdhit]}/cd-hit-est
-i transcripts.fasta
-o transcripts_cdhit.fasta""",
stdfiles[3]:"""{reqd[salmon]}/bin/salmon
index
--index Data/Salmon_Index 
--transcripts transcripts_cdhit.fasta
-p 48""",
stdfiles[4]:"""{reqd[salmon]}/bin/salmon quant
--index Data/Salmon_Index
--libType A
-1 {x}
-2 {y}""",
stdfiles[5]:"""gunzip
-k RNA_SALMON_{ID}/aux_info/eq_classes.txt.gz""",
stdfiles[6]:"""{reqd[hisat2]}/hisat2-build 
transcripts_Corset.fasta
Index""",
stdfiles[7]:"""{reqd[hisat2]}/hisat2
-p {threadsperrun}
--dta

-q
--index Index
-1 {x}
-2 {y}
{ID}.sam"""
}


# Script, Dependencies, Memory 

if __name__ == '__main__':
    app.run(debug=True)
    