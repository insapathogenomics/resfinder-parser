import os
import pytest
import pandas as pd
from resfinder_parser import ResfinderParser, ResfinderCollector


@pytest.mark.parametrize("isolate", ["Isolate01", "Isolate02", "Isolate03", "Isolate04"])
def test_resfinder_parser_basic(test_data_dir, isolate):
    parser = ResfinderParser(test_data_dir, isolate)
    if isolate == "Isolate04":
        assert not parser.has_data
        assert parser.passport.result_summary == "No results found."
        assert parser.resfinder_results_filepath == ""
    else:
        assert parser.has_data
        assert parser.passport.result_summary != ""
        assert os.path.exists(parser.resfinder_results_filepath)
        df = parser.collect_phenotype_results()
        assert isinstance(df, pd.DataFrame)
        assert "isolate_id" in df.columns


@pytest.mark.parametrize("isolate,expected", [
    ("Isolate01", False),
    ("Isolate02", False),
    ("Isolate03", True),
    ("Isolate04", False),
])
def test_resfinder_parser_pointfinder(test_data_dir, isolate, expected):
    parser = ResfinderParser(test_data_dir, isolate)
    assert parser.has_pointfinder_data == expected


def test_resfinder_collector_init(test_data_dir):
    collector = ResfinderCollector(test_data_dir)
    assert set(collector.isolate_dirs) == {"Isolate01", "Isolate02", "Isolate03", "Isolate04"}


def test_resfinder_collector_all_results(test_data_dir):
    collector = ResfinderCollector(test_data_dir)
    results = collector.collect_all_results()
    assert isinstance(results, tuple)
    assert len(results) == 4
    assert isinstance(results[1], pd.DataFrame)
    assert isinstance(results[2], pd.DataFrame)
    assert isinstance(results[3], pd.DataFrame)
