import sys
if sys.version_info < (3, 12):
    raise SystemExit("HPC-T-Assembly requires Python 3.12 or newer.")

from os import getcwd as gc
from os.path import abspath, expandvars, basename, exists, isdir
from os import system
import yaml
from sys import argv
from time import time
import csv
import glob
import io
import json
import os
import re
import shlex
import shutil
import uuid

reqd = {}
base = ['#SBATCH -N ', '#SBATCH -n ', '#SBATCH --mem=', '#SBATCH --account ', '#SBATCH --time ']
Multispecie = False
ExecuteNow = True if argv[1:2] in [["multispecie"], ["Execute"]] else False  # False if argv[1:2] == [] else True


def read_read_pairs(contents):
    """Parse a CSV manifest, ignoring blank lines and validating paired reads."""
    pairs = []
    for row in csv.reader(io.StringIO(contents)):
        if not row or all(not item.strip() for item in row):
            continue
        if len(row) != 2 or not all(item.strip() for item in row):
            raise ValueError("Each manifest row must contain exactly two read paths")
        pairs.append(tuple(abspath(expandvars(item.strip())) for item in row))
    return pairs


def write_read_pairs(path, pairs):
    """Write paired paths as CSV without losing commas or spaces in paths."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerows(pairs)


def discover_read_pairs(directory="Data"):
    """Find paired .fastq/.fq files, including gzip-compressed reads."""
    pattern = re.compile(r"^(.*)_([12])\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)
    found = {}
    for path in sorted(glob.glob(os.path.join(directory, "*"))):
        match = pattern.match(os.path.basename(path))
        if match:
            found.setdefault(match.group(1), {})[match.group(2)] = abspath(path)
    pairs = []
    for sample, mates in found.items():
        if set(mates) != {"1", "2"}:
            raise ValueError(f"Missing mate for read pair {sample!r}")
        pairs.append((mates["1"], mates["2"]))
    if not pairs:
        raise ValueError(f"No paired FASTQ reads found in {directory!r}")
    return pairs


def read_stem(path):
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if path.lower().endswith(suffix):
            return path[:-len(suffix)]
    return os.path.splitext(path)[0]


def move_input_reads(manifest_path="HPC_T_Assembly_Data.txt", target="Data2"):
    """Move every manifest read into target without shell word splitting."""
    os.makedirs(target, exist_ok=True)
    with open(manifest_path, newline="") as f:
        pairs = read_read_pairs(f.read())
    for pair in pairs:
        for source in pair:
            if not exists(source):
                raise FileNotFoundError(f"Input read listed in manifest does not exist: {source}")
            shutil.move(source, target)


def mark_species_complete(species_name, run_id=None):
    """Remove shared software only after every species cleanup succeeds."""
    if basename(species_name) != species_name:
        raise ValueError("Invalid species folder name")
    root = os.path.abspath("..")
    current_context = read_species_cleanup_context()
    context = read_species_cleanup_context(run_id)
    expected = set(context["species"])
    if species_name not in expected:
        raise ValueError(f"Unexpected species cleanup marker: {species_name!r}")
    marker_dir = os.path.join(root, ".species_cleanup", context["run_id"])
    os.makedirs(marker_dir, exist_ok=True)
    marker = os.path.join(marker_dir, species_name + ".done")
    temporary_marker = marker + "." + uuid.uuid4().hex + ".tmp"
    with open(temporary_marker, "w") as f:
        f.write(context["run_id"] + "\n")
    os.replace(temporary_marker, marker)
    completed = {
        name[:-5] for name in os.listdir(marker_dir)
        if name.endswith(".done")
    }
    if expected.issubset(completed) and current_context["run_id"] == context["run_id"]:
        software = os.path.join(root, "Software")
        retired_software = software + ".cleanup-" + uuid.uuid4().hex
        try:
            os.rename(software, retired_software)
        except FileNotFoundError:
            return
        if os.path.islink(retired_software):
            os.unlink(retired_software)
        else:
            shutil.rmtree(retired_software)


def read_species_cleanup_context(run_id=None):
    """Read and validate this species' current shared-software cleanup run."""
    with open(".species_cleanup.json") as f:
        pointer = json.load(f)
    if run_id is not None and pointer.get("run_id") != run_id:
        # Cleanup jobs keep their submitted run ID even if a later run updates
        # this species directory's pointer.
        pointer = {"run_id": run_id}
    run_id = pointer.get("run_id")
    if not isinstance(run_id, str) or not re.fullmatch(r"[a-f0-9]{32}", run_id):
        raise ValueError("Invalid multi-species cleanup metadata; run the top-level setup again")
    manifest_path = os.path.join("..", ".species_cleanup", run_id, "manifest.json")
    with open(manifest_path) as f:
        context = json.load(f)
    if (
        not isinstance(context, dict)
        or context.get("run_id") != run_id
        or not isinstance(context.get("species"), list)
        or not context["species"]
        or any(not isinstance(name, str) or not name or basename(name) != name for name in context["species"])
        or len(set(context["species"])) != len(context["species"])
    ):
        raise ValueError("Invalid multi-species cleanup metadata; run the top-level setup again")
    return context


def write_species_cleanup_context(root, species_names, run_id=None):
    """Write a run-scoped cleanup manifest into each species directory."""
    species_names = list(species_names)
    if not species_names or any(not name or basename(name) != name for name in species_names):
        raise ValueError("Invalid expected species list for cleanup")
    if len(set(species_names)) != len(species_names):
        raise ValueError("Species names must be unique for cleanup")
    run_id = run_id or uuid.uuid4().hex
    if not re.fullmatch(r"[a-f0-9]{32}", run_id):
        raise ValueError("Invalid cleanup run ID")
    context = {"run_id": run_id, "species": species_names}
    marker_dir = os.path.join(root, ".species_cleanup", run_id)
    os.makedirs(marker_dir, exist_ok=True)
    manifest_path = os.path.join(marker_dir, "manifest.json")
    temporary_manifest = manifest_path + "." + uuid.uuid4().hex + ".tmp"
    with open(temporary_manifest, "w") as f:
        json.dump(context, f)
        f.write("\n")
    os.replace(temporary_manifest, manifest_path)
    for species_name in species_names:
        context_path = os.path.join(root, species_name, ".species_cleanup.json")
        temporary_path = context_path + "." + uuid.uuid4().hex + ".tmp"
        with open(temporary_path, "w") as f:
            json.dump(context, f)
            f.write("\n")
        os.replace(temporary_path, context_path)
    return run_id


def clear_species_complete(species_name, run_id=None):
    """Invalidate this species' marker before submitting a new pipeline run."""
    if basename(species_name) != species_name:
        raise ValueError("Invalid species folder name")
    context = read_species_cleanup_context(run_id)
    if species_name not in context["species"]:
        raise ValueError(f"Unexpected species cleanup marker: {species_name!r}")
    marker = os.path.join("..", ".species_cleanup", context["run_id"], species_name + ".done")
    try:
        os.unlink(marker)
    except FileNotFoundError:
        pass


def has_shared_software_context():
    """Whether this species uses shared software owned by its parent directory."""
    return (
        not isdir("Software")
        and isdir("../Software")
        and (exists(".species_cleanup.json") or exists("../.species_cleanup.expected"))
    )


def read_job_specs(config_text):
    specs = []
    for line in config_text.splitlines()[1:]:
        fields = line.split()
        if len(fields) >= 2:
            specs.append((fields[0], fields[1:-1], fields[-1]))
    return specs


def build_submission_script(config_text, start_at=None, shared_software=False):
    """Build submissions, resolving only dependencies included in this run."""
    specs = read_job_specs(config_text)
    if shared_software:
        specs = [spec for spec in specs if spec[0] != "remove_software.sh"]
    if start_at is not None:
        starts = [i for i, spec in enumerate(specs) if spec[0] == start_at]
        if not starts:
            raise ValueError(f"Unknown retry stage: {start_at}")
        specs = specs[starts[0]:]
    active_names = {basename(script).split(".")[0] for script, _, _ in specs}
    commands = []
    for script, dependencies, memory in specs:
        deps = [basename(dep).split(".")[0] for dep in dependencies
                if basename(dep).split(".")[0] in active_names]
        command = f'{basename(script).split(".")[0]}=$(sbatch --parsable'
        if deps:
            deprefs = ["$" + "{" + dep + "}" for dep in deps]
            command += " --dependency=afterany:" + ":".join(deprefs)
        command += f" --mem={memory} {script})"
        commands.append(command)
    return "\n".join(commands) + "\n"


def write_submission_script(start_at=None):
    with open("Config/sbatch.config.txt") as f:
        config_text = f.read()
    shared_software = has_shared_software_context()
    script = build_submission_script(config_text, start_at, shared_software)
    with open("HPC_T_Assembly_Single.sh", "w") as f:
        if shared_software and remove():
            context = read_species_cleanup_context()
            f.write(
                'python HPC_T_Assembly.py clear-species-complete "$(basename "$PWD")" '
                + shlex.quote(context["run_id"]) + " || exit 1\n"
            )
        f.write(script)
        if any(line.startswith("cleanup=") for line in script.splitlines()):
            f.write('printf \'%s\\n\' "$cleanup" > cleanup.jobid\n')
    return shared_software


def submit_jobs(start_at=None):
    write_submission_script(start_at)
    return system("bash HPC_T_Assembly_Single.sh")


def get_config_threads(configf):
    with open(configf) as f:
        for line in f:
            if line.startswith("Threads:"):
                value = int(line.split(":", 1)[1].strip())
                if value < 1:
                    raise ValueError("Configured thread count must be positive")
                return value
    raise ValueError(f"Missing Threads setting in {configf}")


def write_parallel_commands(file_handle, commands, batch_size):
    """Run bounded batches concurrently and wait for every child process."""
    batch_size = max(1, batch_size)
    for start in range(0, len(commands), batch_size):
        for command in commands[start:start + batch_size]:
            file_handle.write(command + " &\n")
        file_handle.write("batch_status=0\n")
        file_handle.write("for child_pid in $(jobs -p); do wait \"$child_pid\" || batch_status=1; done\n")
        file_handle.write("if [ \"$batch_status\" -ne 0 ]; then exit 1; fi\n")


def mainhpc(threads=None):
    reqd = getreqs()  # Get complete path to required software
    with open("HPC_T_Assembly_Data.txt") as f:
        allreads = f.read()
        l = allreads.split("#")
    if len(l) > 1:  # Determine if multiple species
        Multispecie = True
        species = []
        groups = [group for group in allreads.split("#")[1:] if group.strip()]
        for group_number, group in enumerate(groups, 1):
            rows = group.splitlines()
            if not rows or not rows[0].strip():
                continue
            requested_name = "_".join(rows[0].split())
            specie_name = re.sub(r"[^A-Za-z0-9_.-]", "_", requested_name).strip("._")
            if not specie_name:
                specie_name = f"species_{group_number}"
            if specie_name in species:
                raise ValueError(f"Species names must be unique after normalization: {specie_name!r}")
            pairs = read_read_pairs("\n".join(rows[1:]))
            if not pairs:
                raise ValueError(f"No paired reads listed for species {requested_name!r}")
            species.append(specie_name)
            os.makedirs(specie_name, exist_ok=True)
            shutil.copytree("Config", f"{specie_name}/Config", dirs_exist_ok=True)
            shutil.copy2("HPC_T_Assembly.py", specie_name)
            write_read_pairs(
                f"{specie_name}/HPC_T_Assembly_Data.txt",
                [(abspath(expandvars(left)), abspath(expandvars(right))) for left, right in pairs],
            )

        if not species:
            raise ValueError("The multi-species manifest contains no species groups")
        with open("HPC_T_Assembly_Multiple.sh", "w") as genscript:
            if remove():
                expected_species = " ".join(shlex.quote(name) for name in species)
                genscript.write(
                    "python HPC_T_Assembly.py begin-species-cleanup-run "
                    + expected_species + " || exit 1\n"
                )
            for specie_name in species:
                genscript.write(f"cd {specie_name}\npython HPC_T_Assembly.py multispecie\ncd ..\n")
        if remove():
            write_species_cleanup_context(os.path.abspath("."), species)
        if ExecuteNow:  # If execute now --> run multispecie
            s("bash HPC_T_Assembly_Multiple.sh")
        exit()
    pairs = read_read_pairs(allreads)
    if not pairs:
        raise ValueError("HPC_T_Assembly_Data.txt contains no paired reads")
    pairs = [(abspath(expandvars(x)), abspath(expandvars(y))) for x, y in pairs]
    write_read_pairs("HPC_T_Assembly_Data.txt", pairs)
    left = [pair[0] for pair in pairs]
    right = [pair[1] for pair in pairs]
    threads = get_config_threads("Config/fastp.config.txt")
    bowtie_threads = get_config_threads("Config/bowtie2.config.txt")
    salmon_threads = get_config_threads("Config/salmon.config.txt")
    threadsperrun = max(1, bowtie_threads // len(left))
    salmon_threadsperrun = max(1, salmon_threads // len(left))
    fastp_threads = min(threads, 16)
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
    with open(f"pipeline.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(sbatchc)
        f.write(f'cd {gc()}\n')
        fastp_commands = []
        for i in range(len(left)):
            lefti = shlex.quote(left[i])
            righti = shlex.quote(right[i])
            leftio = shlex.quote(read_stem(left[i]))
            rightio = shlex.quote(read_stem(right[i]))
            fastp_commands.append(command.format(**locals()))
        write_parallel_commands(f, fastp_commands, max(1, threads // fastp_threads))

    # Assembly
    lreads = []
    rreads = []

    for i in range(len(left)):
        lreads.append(f"{read_stem(left[i])}_cleaned.fastq")
        rreads.append(f"{read_stem(right[i])}_cleaned.fastq")

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

    # Second Assembler TrinityRnaSeq
    tleft = shlex.quote(",".join(lreads))
    tright = shlex.quote(",".join(rreads))
    trinity = getcommand("Config/trinity.config.txt").format(**locals())
    with open("trinity.sh","w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/trinity.config.txt"))
        f.write(f'cd {gc()}\n')
        f.write(f'\nexport PATH=$PATH:$(pwd)/{reqd["jellyfish"]}\nexport PATH=$PATH:$(pwd)/{reqd["jellyfish"]}/bin\nexport PATH=$PATH:$(pwd)/{reqd["samtools"]}\nexport PATH=$PATH:$(pwd)/{reqd["samtools"]}/bin\nexport PATH=$PATH:$(pwd)/{reqd["bowtie"]}\nexport PATH=$PATH:$(pwd)/{reqd["bowtie"]}/bin\nexport PATH=$PATH:$(pwd)/{reqd["salmon"]}/bin\n')
        # f.write('\nexport PATH=$PATH:$(pwd)/Software/jellyfish-2.3.1\nexport PATH=$PATH:$(pwd)/Software/jellyfish-2.3.1/bin\nexport PATH=$PATH:$(pwd)/Software/samtools-1.21\nexport PATH=$PATH:$(pwd)/Software/samtools-1.21/bin\nexport PATH=$PATH:$(pwd)/Software/bowtie2-2.5.4\nexport PATH=$PATH:$(pwd)/Software/bowtie2-2.5.4/bin\nexport PATH=$PATH:$(pwd)/Software/salmon-latest_linux_x86_64/bin\n')
        f.write(trinity)
        """
cd /g100_scratch/userexternal/tposemar/TestTrinity
export PATH=$PATH:$(pwd)/Software/jellyfish-2.3.1
export PATH=$PATH:$(pwd)/Software/jellyfish-2.3.1/bin
export PATH=$PATH:$(pwd)/Software/samtools-1.21
export PATH=$PATH:$(pwd)/Software/samtools-1.21/bin
export PATH=$PATH:$(pwd)/Software/bowtie2-2.5.4
export PATH=$PATH:$(pwd)/Software/bowtie2-2.5.4/bin
time trinityrnaseq-v2.15.2/Trinity --seqType fq --max_memory 300G --left SRR5759448_1_cleaned.fastq --right SRR5759448_2_cleaned.fastq --CPU 48 --output trinity_out"""

    #Selector

    with open("selector.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write("""#SBATCH -N 1
#SBATCH -n 8
#SBATCH """ + getpartition() + """
#SBATCH --mem=12GB
#SBATCH --time 00:15:00
""" + f"#SBATCH --account {getaccount('Config/assembly.config.txt')}\n" + "#SBATCH -o Selector.out\n")
        f.write(f'cd {gc()}\n')
        selector = r"""perl TRINITY_STATS_PATH/util/TrinityStats.pl trinity_out.Trinity.fasta > trinity_stats.txt
perl TRINITY_STATS_PATH/util/TrinityStats.pl transcripts.fasta > spades_stats.txt

# Extract first N50 values using robust pattern matching
n50_trinity=$(grep -m1 'Contig N50:[[:space:]]*[0-9]\+' trinity_stats.txt | awk '{print $NF}')
n50_spades=$(grep -m1 'Contig N50:[[:space:]]*[0-9]\+' spades_stats.txt | awk '{print $NF}')

# Compare values and output result
if [ "$n50_trinity" -gt "$n50_spades" ]; then
    echo "trinity_stats.txt has higher N50 ($n50_trinity vs $n50_spades)"
    mv transcripts.fasta SpadesTranscripts.fasta
    cp trinity_out.Trinity.fasta transcripts.fasta
elif [ "$n50_spades" -gt "$n50_trinity" ]; then
    echo "spades_stats.txt has higher N50 ($n50_spades vs $n50_trinity)"
else
    echo "Both files have identical N50: $n50_trinity"
fi
""".replace("TRINITY_STATS_PATH", reqd["trinityrnaseq"])
        f.write(selector)

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
        ID = shlex.quote(os.path.basename(x.split("_cleaned.fastq")[0]))
        x = shlex.quote(x)
        y = shlex.quote(y)
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
        write_parallel_commands(f, salmon, max(1, salmon_threads // salmon_threadsperrun))

    with open("salmonpos.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(getsbatch("Config/salmonpos.config.txt"))
        f.write(f'cd {gc()}\n')
        write_parallel_commands(f, salmonpos, get_config_threads("Config/salmonpos.config.txt"))

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
        sample_id = os.path.basename(x.split("_cleaned.fastq")[0])
        sample_id = re.sub(r"_[12]$", "", sample_id)
        ID = shlex.quote(sample_id)
        x = shlex.quote(x)
        y = shlex.quote(y)

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
        write_parallel_commands(f, bowtie, max(1, bowtie_threads // max(1, threadsperrun)))

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
        submit_jobs()


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
    jellyfish
    samtools
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
    retry_field = sbatchc.splitlines()[0].split(",")[-1].strip()
    max_retries = int(retry_field) if retry_field.isdigit() else 0
    multi_species_shared = has_shared_software_context()

    with open("cleanup.sh", "w") as f:
        f.write("""#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH """ + getpartition() + """
#SBATCH --mem=12GB
#SBATCH --time 00:15:00
""" + f"#SBATCH --account {getaccount('Config/assembly.config.txt')}\n" + "#SBATCH -o Verification_Cleaning.out\n")
        if remove() and multi_species_shared:
            context = read_species_cleanup_context()
            f.write(
                'python HPC_T_Assembly.py clear-species-complete "$(basename "$PWD")" '
                + shlex.quote(context["run_id"]) + " || exit 1\n"
            )
        f.write(f"max_retries={max_retries}\n")
        f.write("retry_state=cleanup.retry.state\n[ -f \"$retry_state\" ] || printf '0\\n' > \"$retry_state\"\n")
        f.write("""fix() {
  failed_script=$1
  attempts=$(cat "$retry_state")
  if [[ $attempts -ge $max_retries ]]; then
    echo "Error: Could not fix the issue after $attempts retries."
    exit 1
  fi
  attempts=$((attempts + 1))
  printf '%s\n' "$attempts" > "$retry_state"
  python HPC_T_Assembly.py retry "$failed_script"
  exit 1
}

# Read number of reads from Processes.txt
reads=$(sed -n '2p' Processes.txt | awk -F'|' '{print $2}')

# Verify fastp
cat fastp.err | grep CANCELLED > fastp.verification
if grep -q "CANCELLED" fastp.verification; then
  echo "FastP Verification Failed"
  fix pipeline.sh
  exit 1
fi               
cat fastp.* | grep "Duplication rate:" > fastp.verification
fastp_lines=$(wc -l < fastp.verification)
if [ "$fastp_lines" -lt "$reads" ]; then
  echo "Fastp Verification Failed"
  fix pipeline.sh
  exit 1
fi

# Verify assembly
cat assembly.* | grep "Assembling finished" > assembly.verification
cat assembly.err | grep CANCELLED >> assembly.verification
if grep -q "CANCELLED" assembly.verification; then
  echo "Assembly Verification Failed"
  fix assembly.sh
  exit 1
fi  
if ! grep -q "Assembling finished." assembly.verification; then
  echo "Assembly Verification Failed"
  fix assembly.sh
  exit 1
fi

# Verify CDHIT
cat cdhit.* | grep "writing new database
writing clustering information
program completed" > cdhit.verification
cat cdhit.err | grep CANCELLED >> cdhit.verification
if grep -q "CANCELLED" cdhit.verification; then
  echo "CDHIT Verification Failed"
  fix cdhit.sh
  exit 1
fi  
if ! grep -q "writing new database
writing clustering information
program completed" cdhit.verification; then
  echo "CDHIT Verification Failed"
  fix cdhit.sh
  exit 1
fi

# Verify salmonidx
cat salmonidx.out | grep "Edges construction time:" > salmonidx.verification
cat salmonidx.err | grep CANCELLED >> salmonidx.verification
if grep -q "CANCELLED" salmonidx.verification; then
  echo "SalmonIDX Verification Failed"
  fix salmonidx.sh
  exit 1
fi  
if ! grep -q "Edges construction time:" salmonidx.verification; then
  echo "Salmon idx Verification Failed"
  fix salmonidx.sh
  exit 1
fi

# Verify salmon
cat salmon.err | grep CANCELLED > salmon.verification
if grep -q "CANCELLED" salmon.verification; then
  echo "Salmon Verification Failed"
  fix salmon.sh
  exit 1
fi  
cat salmon.* | grep "done writing equivalence class counts." > salmon.verification
salmon_lines=$(wc -l < salmon.verification)
if [ "$salmon_lines" -lt "$reads" ]; then
  echo "Salmon Verification Failed"
  fix salmon.sh
  exit 1
fi

# Verify corset
cat Corset.* | grep "Finished" > corset.verification
cat Corset.err | grep CANCELLED >> corset.verification
if grep -q "CANCELLED" corset.verification; then
  echo "Corset Verification Failed"
  fix corset.sh
  exit 1
fi  
if ! grep -q "Finished" corset.verification; then
  echo "Corset Verification Failed"
  fix corset.sh
  exit 1
fi
# Verify SalmonPos
cat salmonpos.err | grep "No such file or directory" > salmonpos.verification
if grep -q "No such file or directory" salmonpos.verification; then
    echo "SalmonPos Verification Failed"
    fix salmonpos.sh
    exit 1
fi
# Verify BowtieIDX
cat bowtie2index.err | grep CANCELLED > bowtie.verification
if grep -q "CANCELLED" bowtie.verification; then
  echo "BowtieIDX Verification Failed"
  fix bowtieindex.sh
  exit 1
fi 
cat bowtie2index.err | grep "Renaming Index" > bowtie.verification
if ! grep -q "Renaming Index" bowtie.verification; then
    echo "Bowtie Index Verification Failed"
    fix bowtieindex.sh
    exit 1
fi
  
# Verify Bowtie2
cat bowtie2.err | grep CANCELLED > bowtie2.verification
if grep -q "CANCELLED" bowtie2.verification; then
  echo "Bowtie2 Verification Failed"
  fix bowtie2.sh
  exit 1
fi  
cat bowtie2.* | grep "overall alignment rate" > bowtie.verification
bowtie_lines=$(wc -l < bowtie.verification)
if [ "$bowtie_lines" -lt "$reads" ]; then
  echo "Bowtie2 Verification Failed"
  fix bowtie2.sh
  exit 1
fi

#Verify Busco
cat busco.err | grep CANCELLED > busco.verification
if grep -q "CANCELLED" busco.verification; then
  echo "Busco Verification Failed"
  fix busco.sh
  exit 1
fi
cat busco.* | grep "BUSCO analysis failed!" > busco.verification
if grep -q "BUSCO analysis failed!" busco.verification; then
    echo "Busco Verification Failed"
    fix busco.sh
    exit 1
fi

if grep -q '^transdecoder.sh ' Config/sbatch.config.txt; then
#Verify Transdecoder 
grep "CANCELLED" transdecoder.err > transdecoder.verification || true
if grep -q "CANCELLED" transdecoder.verification
then
  echo "Transdecoder Verification Failed"
  fix transdecoder.sh
  exit 1
fi
cat transdecoder.* | grep "Done preparing long ORFs" > transdecoder.verification
if ! grep -q "Done preparing long ORFs" transdecoder.verification; then
    echo "Transdecoder Verification Failed"
    fix transdecoder.sh
    exit 1
fi

#Verify Transdecoder Predict
grep "transdecoder is finished." transdecoder_predict.err > transdecoder_predict.verification || true
if ! grep -q "transdecoder is finished." transdecoder_predict.verification
then
  echo "Transdecoder Predict Verification Failed"
  fix transdecoder_predict.sh
  exit 1
fi
grep CANCELLED transdecoder_predict.err >> transdecoder_predict.verification || true
if grep -q "CANCELLED" transdecoder_predict.verification
then
  echo "Transdecoder Predict Verification Failed"
  fix transdecoder_predict.sh
  exit 1
fi
fi
""")
        if remove() and not multi_species_shared:
            f.write("bash remove_software.sh || exit 1\n")
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


mkdir -p Data2
python HPC_T_Assembly.py move-inputs || exit 1


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
mv Intermediate_Files/HPC_T_Assembly.py .
mkdir Statistics
mv Intermediate_Files/*stats.txt Statistics 
""")
        if remove() and multi_species_shared:
            context = read_species_cleanup_context()
            f.write(
                'python HPC_T_Assembly.py mark-species-complete "$(basename "$PWD")" '
                + shlex.quote(context["run_id"]) + " || exit 1\n"
            )



def install_missing(name=None, instlist=None):  # Install required software

    installation_instructions = {
        'trinityrnaseq': 'wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.2/trinityrnaseq-v2.15.2.FULL.tar.gz\ntar -zxvf trinityrnaseq-v2.15.2.FULL.tar.gz\ncd trinityrnaseq-v2.15.2\nmake\nmake plugins\ncd ..\n',
        'jellyfin': '\nwget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz\ntar -zxvf jellyfish-2.3.1.tar.gz\ncd jellyfish-2.3.1\n./configure\nmake\ncd ..\n',
        'samtools': 'wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2\ntar -xf samtools-1.21.tar.bz2\ncd samtools-1.21 \n./configure\nmake \ncd ..\n',

        #'trinityrnaseq': '\ngit clone https://github.com/trinityrnaseq/trinityrnaseq.git\ncd trinityrnaseq\nmake\nmake plugins\ncd ..\n',
        'bowtie': 'wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-sra-linux-x86_64.zip -O bowtiebin.zip\nunzip bowtiebin.zip\nmv bowtie2-2.5.4-sra-linux-x86_64 bowtie2-2.5.4\n',
        # 'bowtie': '\nwget https://sourceforge.net/projects/bowtie-bio/files/latest/download -O bowtie.zip\nunzip bowtie.zip\nmkdir bowtie\nmv bowtie*/* bowtie\nrm -rf bowtie2-*\ncd bowtie\ncmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .\n',
        'salmon': '\nwget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gzwget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz\ntar -zxvf salmon-1.10.0_linux_x86_64.tar.gz\n',
        #'samtools': '\ngit clone https://github.com/samtools/samtools.git\ncd samtools\n./configure\nmake\nmake install\ncd ..\n',
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
import yaml
from time import sleep

if __name__ == "__main__":
    if argv[1:2] == ["retry"]:
        if len(argv) < 3:
            raise SystemExit("retry mode requires the failed job script name")
        submit_jobs(argv[2])
        raise SystemExit
    if argv[1:2] == ["move-inputs"]:
        move_input_reads()
        raise SystemExit
    if argv[1:2] == ["mark-species-complete"]:
        if len(argv) < 3:
            raise SystemExit("mark-species-complete requires a species folder name")
        mark_species_complete(argv[2], argv[3] if len(argv) > 3 else None)
        raise SystemExit
    if argv[1:2] == ["clear-species-complete"]:
        if len(argv) < 3:
            raise SystemExit("clear-species-complete requires a species folder name")
        clear_species_complete(argv[2], argv[3] if len(argv) > 3 else None)
        raise SystemExit
    if argv[1:2] == ["begin-species-cleanup-run"]:
        if len(argv) < 3:
            raise SystemExit("begin-species-cleanup-run requires the expected species names")
        write_species_cleanup_context(os.path.abspath("."), argv[2:])
        raise SystemExit
    if os.path.exists("cleanup.retry.state"):
        os.remove("cleanup.retry.state")
    if argv[1:2] == []:
        print("=== Main Menu ===")
        print("1. Execute now\n2. Edit Configuration Batch and Execute Manually")
        if input(": ") == "1":
            ExecuteNow = True
    if argv[1:2] == []:
        install_missing()
    if "HPC_T_Assembly_Data.txt" not in ls():
        write_read_pairs("HPC_T_Assembly_Data.txt", discover_read_pairs())

    write_submission_script()
    with open("HPC_T_Assembly_Data.txt") as f:
        manifest_text = f.read()
    lenleft = 0 if "#" in manifest_text else len(read_read_pairs(manifest_text))

    with open("Processes.txt", "w") as f:
        f.write("Script | Number of Processes\n")
        f.write(
            f"pipeline.sh | {lenleft}\nassembly.sh | 1\nsalmonidx.sh | 1\nsalmon.sh | {lenleft}\nsalmonpos.sh | {lenleft}\ncorset.sh | 1\ncorset2transcript.sh | 1\nbowtieindex.sh | 1\nbowtie2.sh | {lenleft}\ntrasdecoder.sh | 1\ntransdecoder_predict.sh | 1\n")
    cleanup()
    mainhpc()
