import os
import subprocess as sp
import shutil
import unittest
import tempfile

class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.reads_fp = ".tests/reads/"

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output([
            "sunbeam",
            "init",
            "--data_fp",
            self.reads_fp,
            self.project_dir
        ])

        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")
        shutil.copytree(".tests/sunbeam_output", self.output_fp)

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
            "run_demic",
            "--directory",
            self.temp_dir,
            "-n",
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_ptr_fp))
        with open(self.all_ptr_fp) as f:
            self.assertEqual(next(f), "\tSample1\tSample2\tSample3\n")
            for val in next(f).split("\t")[1:]:
                self.assertEqual(round(float(val)), 2)
            for val in next(f).split("\t")[1:]:
                self.assertEqual(round(float(val)), 2)
