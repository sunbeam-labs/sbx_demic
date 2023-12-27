import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile
from pathlib import Path
from typing import Any, Generator


@pytest.fixture
def setup_multi_genome(tmp_path: Path) -> Generator[tuple[Path, Path], Any, None]:
    reads_fp = Path(".tests/data/multi-reads/")

    project_dir = tmp_path / "project"

    sp.check_output(["sunbeam", "init", "--data_fp", str(reads_fp), str(project_dir)])

    yield tmp_path, project_dir


@pytest.fixture
def run_sunbeam_multi_genome(
    setup_multi_genome: tuple[Path, Path]
) -> Generator[tuple[Path, Path], Any, None]:
    tmp_path, project_dir = setup_multi_genome

    output_fp = project_dir / "sunbeam_output"
    logs_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                str(project_dir),
                "all_demic",
                "--directory",
                str(tmp_path),
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(logs_fp, "multi_logs/")
        shutil.copytree(stats_fp, "multi_stats/")
        raise e

    try:
        shutil.copytree(logs_fp, "multi_logs/")
        shutil.copytree(stats_fp, "multi_stats/")
    except FileExistsError:
        pass

    yield output_fp, stats_fp


def test_full_run(run_sunbeam_multi_genome: tuple[Path, Path]) -> None:
    output_fp, stats_fp = run_sunbeam_multi_genome

    all_PTR_fp = output_fp / "mapping" / "demic" / "all_PTR.txt"

    assert all_PTR_fp.exists()
    assert all_PTR_fp.stat().st_size > 0

    with open(all_PTR_fp) as f:
        f.readline()  # Is header
        results = [line.split("\t") for line in f.readlines()]
        print(results)
        assert round(float(results[0][2])) == 2
        assert round(float(results[1][2])) == 3
        assert [float(r[2]) for r in results] == sorted([float(r[2]) for r in results])
