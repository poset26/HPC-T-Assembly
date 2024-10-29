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
    sbatch = ['pipeline.sh  ', 'assembly.sh pipeline.sh ', 'cdhit.sh assembly.sh ', 'salmonidx.sh cdhit.sh ', 'salmon.sh salmonidx.sh ', 'salmonpos.sh salmon.sh ', 'corset.sh salmonpos.sh ', 'corset2transcript.sh corset.sh ', 'statistics.sh corset2transcript.sh ', 'bowtieindex.sh corset2transcript.sh ', 'bowtie2.sh bowtieindex.sh ', 'cleanup.sh bowtie2.sh 12g']
    sbatcht = sbatch[:-1] + ['transdecoder.sh bowtie2.sh ', 'transdecoder_predict.sh transdecoder.sh ', 'cleanup.sh transdecoder_predict.sh 12g']
    stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "bowtie2index", "bowtie2"] #Use amount of threads specified by the user
    singlecore = ["Corset", "Corset2transcript","trinitystats","transdecoder", "transdecoder_predict"] #Use only 1 thread
    addconfigs = ["fastp","assembly","cdhit","Corset","bowtie2","salmon"]  
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
    filelist = [x + ".config.txt" for x in stdfiles + singlecore] + ["sbatch.config.txt"]#, "HPC_T_Assembly.py"]
    #with open("HPC_T_Assembly.py","wb") as f:
        #f.write(hpcpipeline)
    import zipfile
 
    with zipfile.ZipFile('Config.zip', 'w') as zipMe:
        for file in filelist:
            zipMe.write(file, compress_type=zipfile.ZIP_DEFLATED)
    
    return send_file(f"{os.getcwd()}/Config.zip", as_attachment=True)

# file configurations
stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "bowtie2index", "bowtie2"]
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
stdfiles[6]:"""{reqd[bowtie]}/bowtie2-build 
transcripts_Corset.fasta
Index""",
stdfiles[7]:"""{reqd[bowtie]}/bowtie2
-p {threadsperrun}
-x Index
-1 {x}
-2 {y}
-S {ID}.sam
--very-sensitive > log_{ID}.txt"""
}


# Script, Dependencies, Memory 

if __name__ == '__main__':
    app.run(debug=True)
    