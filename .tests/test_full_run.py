import os
import subprocess as sp
import shutil
import unittest
import tempfile

class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.data_fp = os.path.join(self.temp_dir, "data/")
        shutil.copytree(".tests/data/", self.data_fp)

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output([
            "sunbeam",
            "init",
            "--data_fp",
            os.path.join(self.data_fp, "reads/"),
            self.project_dir
        ])

        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")
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
            "all_demic",
            "--directory",
            self.temp_dir,
        ])

        # Check output
        self.assertTrue(os.path.exists(self.all_ptr_fp))
        with open(self.all_ptr_fp) as f:
            self.assertEqual(next(f), "\tTEST0\tTEST1\tTEST2\tTEST3\tTEST4\n")
            for val in next(f).split("\t")[1:]:
                self.assertEqual(round(float(val)), 3)
