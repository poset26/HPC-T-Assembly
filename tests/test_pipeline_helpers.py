import io
import json
import os
import subprocess
import tempfile
import threading
import time
import unittest
from contextlib import contextmanager
from pathlib import Path

import HPC_T_Assembly as pipeline
from HPC_T_Assembly_Config_Utils import busco_config_for_lineage


SBATCH_CONFIG = """# Script, Dependencies, Memory, Retries
pipeline.sh 12g
assembly.sh pipeline.sh 12g
statistics.sh corset2transcript.sh 12g
bowtieindex.sh corset2transcript.sh 12g
busco.sh bowtieindex.sh 12g
bowtie2.sh bowtieindex.sh 12g
transdecoder.sh bowtie2.sh 12g
transdecoder_predict.sh transdecoder.sh 12g
cleanup.sh transdecoder_predict.sh busco.sh statistics.sh 12g
remove_software.sh cleanup.sh 12g
"""


class ManifestTests(unittest.TestCase):
    def test_manifest_accepts_trailing_blank_line_and_quoted_paths(self):
        rows = '\"/reads/sample, one_R1.fastq\",\"/reads/sample, one_R2.fastq\"' + "\n\n"
        self.assertEqual(
            pipeline.read_read_pairs(rows),
            [("/reads/sample, one_R1.fastq", "/reads/sample, one_R2.fastq")],
        )

    def test_discovers_paired_compressed_and_uncompressed_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("sample_1.fastq.gz", "sample_2.fastq.gz", "another_1.fq", "another_2.fq"):
                (root / name).write_text("")
            pairs = pipeline.discover_read_pairs(str(root))
            self.assertEqual(len(pairs), 2)
            self.assertEqual(pipeline.read_stem(str(root / "sample_1.fastq.gz")), str(root / "sample_1"))

    def test_moves_every_pair_without_shell_word_splitting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "reads with spaces"
            source.mkdir()
            names = [
                ("sample, one_R1.fastq", "sample, one_R2.fastq"),
                ("sample two_R1.fastq", "sample two_R2.fastq"),
            ]
            pairs = []
            for pair in names:
                full_pair = []
                for name in pair:
                    path = source / name
                    path.write_text("reads")
                    full_pair.append(str(path))
                pairs.append(tuple(full_pair))
            manifest = root / "manifest.csv"
            pipeline.write_read_pairs(manifest, pairs)
            destination = root / "Data2"
            pipeline.move_input_reads(str(manifest), str(destination))
            self.assertEqual(
                sorted(path.name for path in destination.iterdir()),
                sorted(name for pair in names for name in pair),
            )

    def test_repeated_busco_config_requests_keep_independent_lineages(self):
        original = "busco -l {buscolineage}"
        first = busco_config_for_lineage(original, "metazoa_odb10")
        second = busco_config_for_lineage(original, "vertebrata_odb10")
        self.assertIn("-l metazoa_odb10", first)
        self.assertIn("-l vertebrata_odb10", second)
        self.assertIn("{buscolineage}", original)


class SubmissionTests(unittest.TestCase):
    def test_retry_starts_by_job_name_and_filters_old_dependencies(self):
        script = pipeline.build_submission_script(SBATCH_CONFIG, "busco.sh")
        self.assertTrue(script.startswith("busco=$(sbatch --parsable --mem=12g busco.sh)"))
        self.assertIn("bowtie2=$(sbatch --parsable --mem=12g bowtie2.sh)", script)
        dep = "$" + "{busco}"
        predict_dep = "$" + "{transdecoder_predict}"
        self.assertIn(
            "cleanup=$(sbatch --parsable --dependency=afterany:" + predict_dep + ":" + dep
            + " --mem=12g cleanup.sh)",
            script,
        )
        self.assertNotIn("--dependency=afterany:", script.splitlines()[0])
        self.assertNotIn("--dependency=afterany:}", script)

    def test_fake_sbatch_captures_valid_initial_and_retry_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            fake_sbatch = fake_bin / "sbatch"
            fake_sbatch.write_text(
                "#!/bin/bash\n"
                "printf '%s\\n' \"$*\" >> \"$SBATCH_LOG\"\n"
                "echo 700\n"
            )
            fake_sbatch.chmod(0o755)
            log = root / "sbatch.log"
            env = dict(os.environ, PATH=f"{fake_bin}:{os.environ['PATH']}", SBATCH_LOG=str(log))

            for start_at in (None, "busco.sh"):
                log.write_text("")
                script = pipeline.build_submission_script(SBATCH_CONFIG, start_at)
                subprocess.run(["bash", "-c", script], cwd=root, env=env, check=True)
                calls = log.read_text().splitlines()
                self.assertTrue(calls)
                self.assertTrue(all("afterany:" not in call or "afterany:}" not in call for call in calls))
                self.assertTrue(all("--dependency=afterany:" not in call or "700" in call for call in calls))

    def test_shared_software_is_not_submitted_as_a_species_job(self):
        script = pipeline.build_submission_script(SBATCH_CONFIG, shared_software=True)
        self.assertNotIn("remove_software=", script)
        self.assertIn("cleanup=", script)

    def test_retry_counter_and_orf_checks_are_generated_from_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "Config").mkdir()
            (root / "Config" / "sbatch.config.txt").write_text(
                "# Script, Dependencies, Memory, 2\npipeline.sh 12g\ncleanup.sh pipeline.sh 12g\n"
            )
            (root / "Config" / "assembly.config.txt").write_text(
                "Nodes: 1\nThreads: 1\nMemory: 12GB\nAccount: test\nTime: 00:15:00\n"
                "#Other Sbatch configs\n-p debug\n-o assembly.out\n-e assembly.err\n# assembly\ncommand\n"
            )
            old_cwd = Path.cwd()
            try:
                os.chdir(root)
                pipeline.cleanup()
            finally:
                os.chdir(old_cwd)
            cleanup = (root / "cleanup.sh").read_text()
            self.assertIn("max_retries=2", cleanup)
            self.assertIn("retry_state=cleanup.retry.state", cleanup)
            self.assertIn('python HPC_T_Assembly.py retry \"$failed_script\"', cleanup)
            self.assertIn("if grep -q '^transdecoder.sh ' Config/sbatch.config.txt; then", cleanup)

    def test_retry_state_survives_separate_cleanup_processes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "Config"
            config.mkdir()
            (config / "sbatch.config.txt").write_text(
                "# Script, Dependencies, Memory, 2\npipeline.sh 12g\ncleanup.sh pipeline.sh 12g\n"
            )
            (config / "assembly.config.txt").write_text(
                "Nodes: 1\nThreads: 1\nMemory: 12GB\nAccount: test\nTime: 00:15:00\n"
                "#Other Sbatch configs\n-p debug\n-o assembly.out\n-e assembly.err\n# assembly\ncommand\n"
            )
            (root / "fastp.err").write_text("CANCELLED\n")
            (root / "Processes.txt").write_text("Script | Number of Processes\npipeline.sh | 1\n")
            fake_bin = root / "bin"
            fake_bin.mkdir()
            fake_python = fake_bin / "python"
            fake_python.write_text("#!/bin/bash\nprintf '%s\\n' \"$*\" >> retry-calls.log\n")
            fake_python.chmod(0o755)
            old_cwd = Path.cwd()
            try:
                os.chdir(root)
                pipeline.cleanup()
                env = dict(os.environ, PATH=f"{fake_bin}:{os.environ['PATH']}")
                for expected in ("1", "2"):
                    result = subprocess.run(
                        ["bash", "cleanup.sh"], env=env, check=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    )
                    self.assertEqual(result.returncode, 1)
                    self.assertEqual(Path("cleanup.retry.state").read_text().strip(), expected)
                result = subprocess.run(
                    ["bash", "cleanup.sh"], env=env, check=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                self.assertEqual(result.returncode, 1)
                self.assertEqual(Path("cleanup.retry.state").read_text().strip(), "2")
                self.assertEqual(len(Path("retry-calls.log").read_text().splitlines()), 2)
            finally:
                os.chdir(old_cwd)


class MultiSpeciesCleanupTests(unittest.TestCase):
    def test_old_run_cleanup_does_not_clear_current_run_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            species = root / "species_a"
            species.mkdir()
            old_run = "b" * 32
            current_run = "a" * 32
            pipeline.write_species_cleanup_context(str(root), ["species_a"], run_id=old_run)
            pipeline.write_species_cleanup_context(str(root), ["species_a"], run_id=current_run)
            current_marker = root / ".species_cleanup" / current_run / "species_a.done"
            current_marker.touch()
            old_cwd = Path.cwd()
            try:
                os.chdir(species)
                pipeline.clear_species_complete("species_a", old_run)
            finally:
                os.chdir(old_cwd)
            self.assertTrue(current_marker.exists())

    def test_newer_run_published_before_retirement_keeps_software(self):
        """A run that starts after the stale pointer read must keep Software."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            software = root / "Software"
            software.mkdir()
            (software / "in-use.txt").write_text("active")
            species_a = root / "species_a"
            species_b = root / "species_b"
            species_a.mkdir()
            species_b.mkdir()
            old_run = "b" * 32
            new_run = "c" * 32
            pipeline.write_species_cleanup_context(
                str(root), ["species_a", "species_b"], run_id=old_run
            )
            (root / ".species_cleanup" / old_run / "species_a.done").write_text(old_run + "\n")
            published = {"done": False}
            real_replace = pipeline.os.replace

            def replace_then_publish(src, dst):
                real_replace(src, dst)
                if not published["done"] and os.path.basename(dst) == "species_b.done":
                    published["done"] = True
                    pipeline.write_species_cleanup_context(
                        str(root), ["species_a", "species_b"], run_id=new_run
                    )

            old_cwd = Path.cwd()
            pipeline.os.replace = replace_then_publish
            try:
                os.chdir(species_b)
                pipeline.mark_species_complete("species_b", old_run)
            finally:
                pipeline.os.replace = real_replace
                os.chdir(old_cwd)
            self.assertTrue(published["done"])
            self.assertTrue(software.exists())
            self.assertEqual((software / "in-use.txt").read_text(), "active")
            pointer = json.loads((species_b / ".species_cleanup.json").read_text())
            self.assertEqual(pointer["run_id"], new_run)

    def test_new_run_blocks_until_retirement_decision_finishes(self):
        """Publication waits while the active-run check and rename are in progress."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            software = root / "Software"
            software.mkdir()
            species_a = root / "species_a"
            species_b = root / "species_b"
            species_a.mkdir()
            species_b.mkdir()
            old_run = "d" * 32
            new_run = "e" * 32
            pipeline.write_species_cleanup_context(
                str(root), ["species_a", "species_b"], run_id=old_run
            )
            (root / ".species_cleanup" / old_run / "species_a.done").write_text(old_run + "\n")
            waiting = threading.Event()
            acquired = threading.Event()
            real_lock = pipeline.species_cleanup_lock
            real_listdir = pipeline.os.listdir
            publisher = {"started": False}

            @contextmanager
            def tracking_lock(lock_root):
                is_publisher = threading.current_thread() is publisher.get("thread")
                if is_publisher:
                    waiting.set()
                with real_lock(lock_root):
                    if is_publisher:
                        acquired.set()
                    yield

            def publish():
                pipeline.write_species_cleanup_context(
                    str(root), ["species_a", "species_b"], run_id=new_run
                )

            publisher["thread"] = threading.Thread(target=publish)

            def listdir(path):
                names = real_listdir(path)
                if os.path.basename(path) == old_run and not publisher["started"]:
                    publisher["started"] = True
                    publisher["thread"].start()
                    self.assertTrue(waiting.wait(5), "new run did not reach the cleanup lock")
                    time.sleep(0.2)
                    publisher["blocked"] = publisher["thread"].is_alive() and not acquired.is_set()
                return names

            pipeline.species_cleanup_lock = tracking_lock
            pipeline.os.listdir = listdir
            old_cwd = Path.cwd()
            try:
                os.chdir(species_b)
                pipeline.mark_species_complete("species_b", old_run)
            finally:
                pipeline.species_cleanup_lock = real_lock
                pipeline.os.listdir = real_listdir
                os.chdir(old_cwd)
                if publisher["thread"].is_alive():
                    publisher["thread"].join(5)
            self.assertTrue(publisher.get("blocked"))
            self.assertFalse(publisher["thread"].is_alive())
            pointer = json.loads((species_b / ".species_cleanup.json").read_text())
            self.assertEqual(pointer["run_id"], new_run)
            self.assertFalse(software.exists())

    def test_markers_from_older_runs_do_not_satisfy_current_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            software = root / "Software"
            software.mkdir()
            species_a = root / "species_a"
            species_b = root / "species_b"
            species_a.mkdir()
            species_b.mkdir()
            old_marker_dir = root / ".species_cleanup" / ("b" * 32)
            old_marker_dir.mkdir(parents=True)
            (old_marker_dir / "species_a.done").touch()
            (old_marker_dir / "species_b.done").touch()
            pipeline.write_species_cleanup_context(
                str(root), ["species_a", "species_b"], run_id="a" * 32
            )
            old_cwd = Path.cwd()
            try:
                os.chdir(species_a)
                pipeline.mark_species_complete("species_a")
                self.assertTrue(software.exists())
                os.chdir(species_b)
                pipeline.mark_species_complete("species_b")
            finally:
                os.chdir(old_cwd)
            self.assertFalse(software.exists())

    def test_shared_software_remains_until_every_species_finishes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            software = root / "Software"
            software.mkdir()
            species_a = root / "species_a"
            species_b = root / "species_b"
            species_a.mkdir()
            species_b.mkdir()
            pipeline.write_species_cleanup_context(
                str(root), ["species_a", "species_b"], run_id="a" * 32
            )
            old_cwd = Path.cwd()
            try:
                os.chdir(species_a)
                pipeline.mark_species_complete("species_a")
                self.assertTrue(software.exists())
                os.chdir(species_b)
                pipeline.mark_species_complete("species_b")
            finally:
                os.chdir(old_cwd)
            self.assertFalse(software.exists())

    def test_shared_software_marker_is_last_cleanup_action(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "Software").mkdir()
            species = root / "species_a"
            (species / "Config").mkdir(parents=True)
            (root / "species_b").mkdir()
            pipeline.write_species_cleanup_context(
                str(root), ["species_a", "species_b"], run_id="a" * 32
            )
            (species / "Config" / "sbatch.config.txt").write_text(
                "# Script, Dependencies, Memory, 0\nremove_software.sh cleanup.sh 12g\n"
            )
            (species / "Config" / "assembly.config.txt").write_text(
                "Nodes: 1\nThreads: 1\nMemory: 12GB\nAccount: test\nTime: 00:15:00\n"
                "#Other Sbatch configs\n-p debug\n-o assembly.out\n-e assembly.err\n# assembly\ncommand\n"
            )
            old_cwd = Path.cwd()
            try:
                os.chdir(species)
                pipeline.cleanup()
            finally:
                os.chdir(old_cwd)
            generated = (species / "cleanup.sh").read_text()
            self.assertIn('clear-species-complete "$(basename "$PWD")" ', generated)
            self.assertGreater(
                generated.index("mark-species-complete"),
                generated.index("mv Intermediate_Files/*stats.txt Statistics"),
            )


class ParallelCommandTests(unittest.TestCase):
    def test_generated_batches_wait_for_every_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp) / "completed"
            script = io.StringIO()
            pipeline.write_parallel_commands(
                script,
                [f"(sleep 0.1; echo first >> {marker})", f"(sleep 0.3; echo second >> {marker})"],
                2,
            )
            script.write(f"test \"$(wc -l < {marker})\" -eq 2\n")
            subprocess.run(["bash", "-c", script.getvalue()], check=True)

    def test_generated_batches_propagate_child_failure(self):
        script = io.StringIO()
        pipeline.write_parallel_commands(script, ["true", "false"], 2)
        result = subprocess.run(["bash", "-c", script.getvalue()])
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
