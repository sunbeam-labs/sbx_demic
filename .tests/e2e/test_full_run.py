import csv
import os
import pytest
import shutil
import subprocess as sp
import tempfile


@pytest.fixture
def setup():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

    config_str = f"sbx_demic: {{extras: '-W 701 -D 50 -L 31'}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    output_fp = os.path.join(project_dir, "sunbeam_output")

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                project_dir,
                "all_demic",
                "--directory",
                temp_dir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
        sp.CalledProcessError(e)

    all_ptr_fp = os.path.join(output_fp, "mapping/demic/DEMIC_OUT/all_PTR.txt")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield all_ptr_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    all_ptr_fp, benchmarks_fp = run_sunbeam

    # Check output
    assert os.path.exists(all_ptr_fp)
    with open(all_ptr_fp) as f:
        assert next(f) == "\tTEST0\tTEST1\tTEST2\tTEST3\tTEST4\n"
        for val in next(f).split("\t")[1:]:
            assert round(float(val)) == 3
