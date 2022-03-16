import os
import sys

import subprocess as sp
import shutil
import unittest
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.data_fp = os.path.join(self.temp_dir, "data/")
        shutil.copytree(".tests/data/", self.data_fp)

        self.samples_fp = os.path.join(self.temp_dir, "samples.csv")
        self.samples_content = """
TEST0,{r0r1},{r0r2}\n
TEST1,{r1r1},{r1r2}\n
TEST2,{r2r1},{r2r2}\n
TEST3,{r3r1},{r3r2}\n
TEST4,{r4r1},{r4r2}\n
""".format(
                r0r1 = os.path.join(self.data_fp, "reads/TEST0_R1.fastq.gz"),
                r0r2 = os.path.join(self.data_fp, "reads/TEST0_R2.fastq.gz"),
                r1r1 = os.path.join(self.data_fp, "reads/TEST1_R1.fastq.gz"),
                r1r2 = os.path.join(self.data_fp, "reads/TEST1_R2.fastq.gz"),
                r2r1 = os.path.join(self.data_fp, "reads/TEST2_R1.fastq.gz"),
                r2r2 = os.path.join(self.data_fp, "reads/TEST2_R2.fastq.gz"),
                r3r1 = os.path.join(self.data_fp, "reads/TEST3_R1.fastq.gz"),
                r3r2 = os.path.join(self.data_fp, "reads/TEST3_R2.fastq.gz"),
                r4r1 = os.path.join(self.data_fp, "reads/TEST4_R1.fastq.gz"),
                r4r2 = os.path.join(self.data_fp, "reads/TEST4_R2.fastq.gz"))
        with open(self.samples_fp, "w") as f:
            f.write(self.samples_content)
        
        self.config_fp = os.path.join(self.temp_dir, "sunbeam_config.yml")
        self.config_content = """
            all:\n
              root: {root}\n
              output_fp: sunbeam_output\n
              samplelist_fp: samples.csv\n
              paired_end: true\n
              download_reads: false\n
              version: 2.1.1+dev0.g43432e1.d20220208\n
            qc:\n
              suffix: qc\n
              threads: 4\n
              java_heapsize: 512M\n
              leading: 3\n
              trailing: 3\n
              slidingwindow: [4, 15]\n
              minlen: 36\n
              adapter_fp: ''\n
              fwd_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
              rev_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
              kz_threshold: 0.55\n
              pct_id: 0.5\n
              frac: 0.6\n
              host_fp: ''\n
            classify:\n
              suffix: classify\n
              threads: 4\n
              kraken_db_fp: ''\n
            assembly:\n
              suffix: assembly\n
              min_length: 300\n
              threads: 4\n
            annotation:\n
              suffix: annotation\n
              min_contig_len: 500\n
              circular_kmin: 10\n
              circular_kmax: 1000\n
              circular_min_len: 3500\n
            blast:\n
              threads: 4\n
            blastdbs:\n
              root_fp: ''\n
            mapping:\n
              suffix: mapping\n
              genomes_fp: ''\n
              samtools_opts: ''\n
              threads: 4\n
            download:\n
              suffix: download\n
              threads: 4\n
              sbx_coassembly: ''\n
              threads: 4\n
              group_file: ''\n
            sbx_coassembly:\n
              threads: 4\n
              group_file: ''\n
            sbx_demic:\n
              threads: 4 #--thread_num\n
              keepall: FALSE #--output_all\n
              extras: "-W 701 -D 50 -L 31" # Other parameters passed to DEMIC.pl\n
            """.format(root = self.temp_dir)
        with open(self.config_fp, "w") as f:
            f.write(self.config_content)

        self.output_fp = os.path.join(self.temp_dir, "sunbeam_output")
        #shutil.copytree(".tests/data/sunbeam_output", self.output_fp)

        self.all_ptr_fp = os.path.join(self.output_fp, "mapping/demic/DEMIC_OUT/all_PTR.txt")
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_full_run(self):
        # Run the test job.
        sp.check_output([
            "sunbeam",
            "run",
             "--configfile", 
            self.config_fp,
            "--use-conda",
            "all_demic",
            "--directory",
            self.temp_dir,
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_ptr_fp))
        with open(self.all_ptr_fp) as f:
            self.assertEqual(next(f), "\tTEST0\tTEST1\tTEST2\tTEST3\tTEST4\n")
            for val in next(f).split("\t")[1:]:
                self.assertEqual(int(val), 3)