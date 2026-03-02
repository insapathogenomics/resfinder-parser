import os
import pytest
import pandas as pd
from resfinder_parser import ResfinderParser, ResfinderCollector


@pytest.mark.parametrize("isolate", ["Isolate01", "Isolate02", "Isolate03", "Isolate04"])
def test_resfinder_parser_basic(test_data_dir, isolate):
    parser = ResfinderParser(test_data_dir, isolate)
    if isolate == "Isolate04":
        assert not parser.has_json_data
        assert parser.passport.result_summary == "No results found."
        assert parser.resfinder_results_filepath == ""
    else:
        assert parser.has_json_data
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
    assert set(collector.isolate_dirs) == {"Isolate01", "Isolate02", "Isolate03", "Isolate04", "Isolate05"}


def test_resfinder_collector_all_results(test_data_dir):
    collector = ResfinderCollector(test_data_dir)
    collector.inspect_directories()
    results = collector.collect_all_results()
    assert isinstance(results, tuple)
    assert len(results) == 4
    assert isinstance(results[1], pd.DataFrame)
    assert isinstance(results[2], pd.DataFrame)
    assert isinstance(results[3], pd.DataFrame)


def test_collect_variations_results_basic(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate03")
    assert parser.has_json_data
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)
    expected_columns = [
        "Mutation", "Nucleotide change", "Amino acid change",
        "Resistance", "PMID", "Databases", "isolate_id", "analysis_date", "version"
    ]
    for col in expected_columns:
        assert col in df.columns
    assert len(df) > 0


def test_collect_variations_results_no_data(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate01")
    assert parser.has_json_data
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_collect_variations_results_no_json(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate04")
    assert not parser.has_json_data
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)


def test_collect_variations_results_with_exclude(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate05", exclude_databases=["pointfinder"])
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)
    assert "PointFinder" not in df["Databases"].values
    assert len(df) == 1
    assert "DisinFinder" in df["Databases"].values[0]


def test_exclude_databases_in_parser(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate05", exclude_databases=["disinfinder"])
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert "DisinFinder" not in "; ".join(df["Databases"].values)


def test_exclude_databases_case_insensitive(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate05", exclude_databases=["POINTFINDER"])
    df = parser.collect_variations_results()
    assert isinstance(df, pd.DataFrame)
    assert "PointFinder" not in df["Databases"].values
    assert len(df) > 0


def test_exclude_databases_in_collector(test_data_dir):
    collector = ResfinderCollector(test_data_dir, exclude_databases=["pointfinder"])
    collector.inspect_directories()
    region_results, variation_results, isolate_summaries, combined = collector.collect_all_results()
    if not variation_results.empty:
        assert "PointFinder" not in variation_results["Databases"].values


def test_collector_with_exclude_and_variations(test_data_dir):
    collector = ResfinderCollector(test_data_dir, exclude_databases=["disinfinder"])
    collector.inspect_directories()
    region_results, variation_results, isolate_summaries, combined = collector.collect_all_results()
    if not variation_results.empty:
        variation_from_disin = variation_results[
            variation_results["Databases"].str.contains("DisinFinder", na=False)
        ]
        assert len(variation_from_disin) == 0


def test_exclude_databases_in_seq_regions(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate05", exclude_databases=["resfinder"])
    seq_reg_df = parser.seq_regions_parse()
    assert isinstance(seq_reg_df, pd.DataFrame)
    assert "ResFinder" not in "; ".join(seq_reg_df["databases"].values)


def test_exclude_databases_in_summary(test_data_dir):
    parser = ResfinderParser(test_data_dir, "Isolate05", exclude_databases=["disinfinder"])
    summary = parser.isolate_summary()
    assert isinstance(summary, pd.DataFrame)
    assert "DisinFinder" not in summary["databases"].values[0]

    