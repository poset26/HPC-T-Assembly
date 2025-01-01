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
    sbatch = ['pipeline.sh  ', 'assembly.sh pipeline.sh ', 'cdhit.sh assembly.sh ', 'salmonidx.sh cdhit.sh ', 'salmon.sh salmonidx.sh ', 'salmonpos.sh salmon.sh ', 'corset.sh salmonpos.sh ', 'corset2transcript.sh corset.sh ', 'statistics.sh corset2transcript.sh ', 'bowtieindex.sh corset2transcript.sh ', 'busco.sh bowtieindex.sh ','bowtie2.sh bowtieindex.sh ', 'cleanup.sh bowtie2.sh 12g']
    sbatcht = sbatch[:-1] + ['transdecoder.sh bowtie2.sh ', 'transdecoder_predict.sh transdecoder.sh ', 'cleanup.sh transdecoder_predict.sh 12g']
    stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "bowtie2index", "bowtie2", "busco"] #Use amount of threads specified by the user
    singlecore = ["Corset", "Corset2transcript","trinitystats","transdecoder", "transdecoder_predict"] #Use only 1 thread
    addconfigs = ["fastp","assembly","cdhit","Corset","bowtie2","salmon"]  
    sbatch = sbatcht if request.form.get("ORFs") == "true" else sbatch #Add ORF prediction to the pipeline if requested
    buscolineage = request.form["busco_additional_configurations"] if request.form["busco_additional_configurations"] in lineages else "nematoda_odb10"
    
    defaultconfigs["busco"] = defaultconfigs["busco"].replace("{buscolineage}",buscolineage)
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
        if request.form.get("retries") != "0":
            f.write(f"# Script, Dependencies, Memory, {request.form.get("retries")}\n")
        else:
            f.write("# Script, Dependencies, Memory, Retries \n")
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
stdfiles = ["fastp", "assembly", "cdhit" ,"salmonidx", "salmon", "salmonpos", "bowtie2index", "bowtie2","busco"]
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
-t {threads}
-o ASSEMBLY/RNA-Spades""",
stdfiles[2]:"""{reqd[cdhit]}/cd-hit-est
-i transcripts.fasta
-o transcripts_cdhit.fasta""",
stdfiles[3]:"""{reqd[salmon]}/bin/salmon
index
--index Data/Salmon_Index 
--transcripts transcripts_cdhit.fasta
-p {threads}""",
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
--very-sensitive > log_{ID}.txt""",
stdfiles[8]:"""{reqd[busco]}/bin/busco 
-i transcripts_Corset.fasta 
-l {buscolineage}
-o busco_nematoda 
-m transcriptome -c {threads}"""
}

lineages = ['bacteria_odb10', 'acidobacteria_odb10', 'actinobacteria_phylum_odb10', 'actinobacteria_class_odb10', 'corynebacteriales_odb10', 'micrococcales_odb10', 'propionibacteriales_odb10', 'streptomycetales_odb10', 'streptosporangiales_odb10', 'coriobacteriia_odb10', 'coriobacteriales_odb10', 'aquificae_odb10', 'bacteroidetes-chlorobi_group_odb10', 'bacteroidetes_odb10', 'bacteroidia_odb10', 'bacteroidales_odb10', 'cytophagia_odb10', 'cytophagales_odb10', 'flavobacteriia_odb10', 'flavobacteriales_odb10', 'sphingobacteriia_odb10', 'chlorobi_odb10', 'chlamydiae_odb10', 'chloroflexi_odb10', 'cyanobacteria_odb10', 'chroococcales_odb10', 'nostocales_odb10', 'oscillatoriales_odb10', 'synechococcales_odb10', 'firmicutes_odb10', 'bacilli_odb10', 'bacillales_odb10', 'lactobacillales_odb10', 'clostridia_odb10', 'clostridiales_odb10', 'thermoanaerobacterales_odb10', 'selenomonadales_odb10', 'tissierellia_odb10', 'tissierellales_odb10', 'fusobacteria_odb10', 'fusobacteriales_odb10', 'planctomycetes_odb10', 'proteobacteria_odb10', 'alphaproteobacteria_odb10', 'rhizobiales_odb10', 'rhizobium-agrobacterium_group_odb10', 'rhodobacterales_odb10', 'rhodospirillales_odb10', 'rickettsiales_odb10', 'sphingomonadales_odb10', 'betaproteobacteria_odb10', 'burkholderiales_odb10', 'neisseriales_odb10', 'nitrosomonadales_odb10', 'delta-epsilon-subdivisions_odb10', 'deltaproteobacteria_odb10', 'desulfobacterales_odb10', 'desulfovibrionales_odb10', 'desulfuromonadales_odb10', 'epsilonproteobacteria_odb10', 'campylobacterales_odb10', 'gammaproteobacteria_odb10', 'alteromonadales_odb10', 'cellvibrionales_odb10', 'chromatiales_odb10', 'enterobacterales_odb10', 'legionellales_odb10', 'oceanospirillales_odb10', 'pasteurellales_odb10', 'pseudomonadales_odb10', 'thiotrichales_odb10', 'vibrionales_odb10', 'xanthomonadales_odb10', 'spirochaetes_odb10', 'spirochaetia_odb10', 'spirochaetales_odb10', 'synergistetes_odb10', 'tenericutes_odb10', 'mollicutes_odb10', 'entomoplasmatales_odb10', 'mycoplasmatales_odb10', 'thermotogae_odb10', 'verrucomicrobia_odb10', 'archaea_odb10', 'thaumarchaeota_odb10', 'thermoprotei_odb10', 'thermoproteales_odb10', 'sulfolobales_odb10', 'desulfurococcales_odb10', 'euryarchaeota_odb10', 'thermoplasmata_odb10', 'methanococcales_odb10', 'methanobacteria_odb10', 'methanomicrobia_odb10', 'methanomicrobiales_odb10', 'halobacteria_odb10', 'halobacteriales_odb10', 'natrialbales_odb10', 'haloferacales_odb10', 'eukaryota_odb10', 'alveolata_odb10', 'apicomplexa_odb10', 'aconoidasida_odb10', 'plasmodium_odb10', 'coccidia_odb10', 'euglenozoa_odb10', 'fungi_odb10', 'ascomycota_odb10', 'dothideomycetes_odb10', 'capnodiales_odb10', 'pleosporales_odb10', 'eurotiomycetes_odb10', 'chaetothyriales_odb10', 'eurotiales_odb10', 'onygenales_odb10', 'leotiomycetes_odb10', 'helotiales_odb10', 'saccharomycetes_odb10', 'sordariomycetes_odb10', 'glomerellales_odb10', 'hypocreales_odb10', 'basidiomycota_odb10', 'agaricomycetes_odb10', 'agaricales_odb10', 'boletales_odb10', 'polyporales_odb10', 'tremellomycetes_odb10', 'microsporidia_odb10', 'mucoromycota_odb10', 'mucorales_odb10', 'metazoa_odb10', 'arthropoda_odb10', 'arachnida_odb10', 'insecta_odb10', 'endopterygota_odb10', 'diptera_odb10', 'hymenoptera_odb10', 'lepidoptera_odb10', 'hemiptera_odb10', 'mollusca_odb10', 'nematoda_odb10', 'vertebrata_odb10', 'actinopterygii_odb10', 'cyprinodontiformes_odb10', 'tetrapoda_odb10', 'mammalia_odb10', 'eutheria_odb10', 'euarchontoglires_odb10', 'glires_odb10', 'primates_odb10', 'laurasiatheria_odb10', 'carnivora_odb10', 'cetartiodactyla_odb10', 'sauropsida_odb10', 'aves_odb10', 'passeriformes_odb10', 'stramenopiles_odb10', 'viridiplantae_odb10', 'chlorophyta_odb10', 'embryophyta_odb10', 'liliopsida_odb10', 'poales_odb10', 'eudicots_odb10', 'brassicales_odb10', 'fabales_odb10', 'solanales_odb10', 'dataset)', 'alphaherpesvirinae_odb10', 'baculoviridae_odb10', 'rudiviridae_odb10', 'betaherpesvirinae_odb10', 'herpesviridae_odb10', 'poxviridae_odb10', 'tevenvirinae_odb10', 'aviadenovirus_odb10', 'enquatrovirus_odb10', 'teseptimavirus_odb10', 'bclasvirinae_odb10', 'fromanvirus_odb10', 'skunavirus_odb10', 'betabaculovirus_odb10', 'pahexavirus_odb10', 'alphabaculovirus_odb10', 'tunavirinae_odb10', 'simplexvirus_odb10', 'gammaherpesvirinae_odb10', 'varicellovirus_odb10', 'cheoctovirus_odb10', 'guernseyvirinae_odb10', 'tequatrovirus_odb10', 'chordopoxvirinae_odb10', 'peduovirus_odb10', 'iridoviridae_odb10', 'spounavirinae_odb10']
# Script, Dependencies, Memory 

if __name__ == '__main__':
    app.run(debug=True)
    