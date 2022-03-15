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
                r0r1 = os.path.join(self.temp_dir, "data/TEST0_R1.fastq.gz"),
                r0r2 = os.path.join(self.temp_dir, "data/TEST0_R2.fastq.gz"),
                r1r1 = os.path.join(self.temp_dir, "data/TEST1_R1.fastq.gz"),
                r1r2 = os.path.join(self.temp_dir, "data/TEST1_R2.fastq.gz"),
                r2r1 = os.path.join(self.temp_dir, "data/TEST2_R1.fastq.gz"),
                r2r2 = os.path.join(self.temp_dir, "data/TEST2_R2.fastq.gz"),
                r3r1 = os.path.join(self.temp_dir, "data/TEST3_R1.fastq.gz"),
                r3r2 = os.path.join(self.temp_dir, "data/TEST3_R2.fastq.gz"),
                r4r1 = os.path.join(self.temp_dir, "data/TEST4_R1.fastq.gz"),
                r4r2 = os.path.join(self.temp_dir, "data/TEST4_R2.fastq.gz"))
        with open(self.samples_fp, "w") as f:
            f.write(self.samples_content)
        
        self.config_fp = os.path.join(self.temp_dir, "sunbeam_config.yml")
        self.config_content = """
            all:\n
            \troot: {root}\n
            \toutput_fp: sunbeam_output\n
            \tsamplelist_fp: samples.csv\n
            \tpaired_end: true\n
            \tdownload_reads: false\n
            \tversion: 2.1.1+dev0.g43432e1.d20220208\n
            qc:\n
            \tsuffix: qc\n
            \tthreads: 4\n
            \tjava_heapsize: 512M\n
            \tleading: 3\n
            \ttrailing: 3\n
            \tslidingwindow: [4, 15]\n
            \tminlen: 36\n
            \tadapter_fp: \n
            \tfwd_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
            \trev_adapters: [GTTTCCCAGTCACGATC, GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC]\n
            \tkz_threshold: 0.55\n
            \tpct_id: 0.5\n
            \tfrac: 0.6\n
            \thost_fp: ''\n
            classify:\n
            \tsuffix: classify\n
            \tthreads: 4\n
            \tkraken_db_fp: ''\n
            assembly:\n
            \tsuffix: assembly\n
            \tmin_length: 300\n
            \tthreads: 4\n
            annotation:\n
            \tsuffix: annotation\n
            \tmin_contig_len: 500\n
            \tcircular_kmin: 10\n
            \tcircular_kmax: 1000\n
            \tcircular_min_len: 3500\n
            blast:\n
            \tthreads: 4\n
            \tblastdbs:\n
            \troot_fp: ''\n
            mapping:\n
            \tsuffix: mapping\n
            \tgenomes_fp: {ref}\n
            \tsamtools_opts: ''\n
            \tthreads: 4\n
            download:\n
            \tsuffix: download\n
            \tthreads: 4\n
            \tsbx_coassembly:\n
            \tthreads: 4\n
            \tgroup_file: ''\n
            sbx_demic:\n
            \tthreads: 4 #--thread_num\n
            \tkeepall: FALSE #--output_all\n
            \textras: "-W 701 -D 50 -L 31" # Other parameters passed to DEMIC.pl\n
            """.format(root = os.path.join(self.temp_dir),
                ref = os.path.join(self.data_fp, "data/reference/akk-genome.fasta"))
        with open(self.config_fp, "w") as f:
            f.write(self.config_content)

        self.output_fp = os.path.join(self.temp_dir, "sunbeam_output")
        shutil.copytree(".tests/data/sunbeam_output", self.output_fp)

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
            "--"
            "--directory",
            self.temp_dir,
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_ptr_fp))
        with open(self.all_ptr_fp) as f:
            self.assertEqual(next(f), "\tTEST0\tTEST1\tTEST2\tTEST3\tTEST4")
            for val in next(f).split("\t")[1:]:
                self.assertEqual(int(val), 10)