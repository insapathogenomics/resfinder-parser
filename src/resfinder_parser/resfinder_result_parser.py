from typing import Optional
import pandas as pd
import json
import os
import datetime
import logging
from typing import List
from .data_classes import IsolatePhenotypes, Phenotype, IsolateSummary


class ResfinderParser:
    pointfinder_results_filename: str = "PointFinder_results.txt"
    resfinder_results_suffix = ".json"
    databases_to_exclude = []

    def __init__(self, RESFINDER_dir, isolate_dir, exclude_databases=None):
        if exclude_databases is None:
            exclude_databases = []
        self.databases_to_exclude = exclude_databases
        self.isolate_id = isolate_dir
        
        if os.path.isdir(os.path.join(RESFINDER_dir, isolate_dir, "resfinder_results")):
            isolate_dir = os.path.join(isolate_dir, "resfinder_results")

        self.pointfinder_results_filepath = os.path.join(
            RESFINDER_dir, isolate_dir, self.pointfinder_results_filename
        )

        json_files = [
            f
            for f in os.listdir(os.path.join(RESFINDER_dir, isolate_dir))
            if f.endswith(self.resfinder_results_suffix)
        ]

        if len(json_files) == 0:
            logging.info(
                f"No resfinder json files found for isolate {self.isolate_id}."
            )
            

        if len(json_files) > 1:
            logging.warning(
                f"Multiple resfinder json files found for isolate {self.isolate_id}. Using {json_files[0]}"
            )
        


        self.time_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        self.has_json_data = len(json_files) > 0

        if self.has_json_data:
            self.resfinder_results_filepath = os.path.join(
                RESFINDER_dir, isolate_dir, json_files[0]
            )
            self.passport = self.collect_summary()

        else:
            self.passport = self.empty_passport
            self.resfinder_results_filepath = ""


        self.has_pointfinder_data = (
            os.path.isfile(self.pointfinder_results_filepath)
            and os.path.getsize(self.pointfinder_results_filepath) > 0
        )

        self.phenotypes = IsolatePhenotypes(
            self.isolate_id, self.passport.result_summary
        )

    @property
    def has_resfinder_data(self):
        return self.has_json_data

    def json_parse_antibiotics(self):
        """Parse the resfinder json output into a pandas dataframe."""

        if self.has_json_data is False:
            return self 
        
        content = self.read_json()
        ## info to get:

        for phenotype_str, data in content["phenotypes"].items():
            antibiotic = Phenotype(phenotype_str)
            antibiotic.feed_data(data)
            self.phenotypes.add_phenotype(antibiotic)

    def resfinder_json_summary(self, content) -> IsolateSummary:
        """
        Parse the resfinder json output into a IsolateSummary object."""
        analysis_key = content["key"]
        provided_species = content["provided_species"]
        result_summary = content["result_summary"]
        databases: dict = content["databases"]
        databases = "; ".join(
            [
                f"{v['database_name']}-{v['database_version']}"
                for k, v in databases.items()
                if v["database_name"].lower() not in self.databases_to_exclude
            ]
        )
        return IsolateSummary(analysis_key, provided_species, result_summary, databases)

    @property
    def empty_passport(self):
        return IsolateSummary("", "", "No results found.", "")

    def read_json(self):
        with open(self.resfinder_results_filepath, "r") as f:
            results = json.load(f)

        return results

    def collect_summary(self) -> IsolateSummary:
        results = self.read_json()
        return self.resfinder_json_summary(results)

    def add_analysis_columns(self, df: pd.DataFrame):
        df["isolate_id"] = self.isolate_id
        df["analysis_date"] = self.time_now
        df["version"] = self.passport.key

        # place the isolate_id, analysis_date and resfinder version columns at the front
        cols = df.columns.tolist()
        cols = cols[-3:] + cols[:-3]
        df = df[cols]

        return df

    def collect_pointfinder_results(self):
        results = pd.read_csv(self.pointfinder_results_filepath, sep="\t")

        results = self.add_analysis_columns(results)

        return results

    def collect_phenotype_results(self):

        if self.has_json_data is False:
            return self.empty_isolate_summary

        self.json_parse_antibiotics()
        isolate_results = self.phenotypes.phenotype_dataframe()
        isolate_results = self.add_analysis_columns(isolate_results)

        return isolate_results
    
    def collect_variations_results(self):
        results = self.read_json()

        seq_var = []
        for var, map_dict in results["seq_variations"].items():
            
            ## filter out if all databases are in the exclude list
            databases = [map_dict["ref_database"]]
            databases_simple = [x.split("-")[0].lower() for x in databases]

            if all(db in self.databases_to_exclude for db in databases_simple):
                continue
            
            phenotypes = "Unknown"
            if map_dict["phenotypes"] is not None and len(map_dict["phenotypes"]) > 0:
                phenotypes = ",".join(map_dict["phenotypes"])
            var_id = var.split(";;")[0]
            protein_change = 'RNA mutations'
            if 'ref_aa' in map_dict:
                protein_change = f"{map_dict['ref_aa'].upper()} -> {map_dict['var_aa'].upper()}"
            
            pmids = "-"
            if 'pmids' in map_dict and map_dict['pmids'] is not None and len(map_dict['pmids']) > 0:
                pmids = ",".join(map_dict['pmids'])

            databases = ",".join(databases)

            for _ in map_dict['seq_regions']:
                seq_var.append([
                    " ".join([var_id, map_dict["seq_var"]]),
                    f"{map_dict['ref_codon'].upper()} -> {map_dict['var_codon'].upper()}",
                    protein_change,
                    phenotypes,
                    pmids,
                    databases
                ])

        seq_var = pd.DataFrame(seq_var, columns=[
            "Mutation",
            "Nucleotide change",
            "Amino acid change",
            "Resistance",
            "PMID",
            "Databases"
        ])

        self.add_analysis_columns(seq_var)

        return seq_var

    def seq_regions_parse(self):
        results = self.read_json()

        seq_reg = []

        for reg, map in results["seq_regions"].items():

            ## filter out if all databases are in the exclude list
            databases = map["ref_database"]
            databases_simple = [x.split("-")[0].lower() for x in databases]

            if all(db in self.databases_to_exclude for db in databases_simple) is True:
                continue

            databases = ",".join(databases)

            for phenotype in map["phenotypes"]:
                seq_reg.append([reg, phenotype, map["query_id"], map["identity"], map['coverage'], map['grade'], databases])

        seq_reg_df = pd.DataFrame(
            seq_reg, columns=["seq_region", "antibiotic", "query_id", "identity", "coverage", "grade", "databases"]
        )

        return seq_reg_df

    def extend_seq_region_results(self, resfinder_df):

        seq_reg_df = self.seq_regions_parse()

        def contig_id_info(phenotype: pd.DataFrame):
            pheno_df = seq_reg_df[seq_reg_df.antibiotic == phenotype]

            if pheno_df.shape[0] == 0:
                return pd.Series(["", 0, 0, 0, ""])

            output = pd.Series(
                ["; ".join(pheno_df.query_id.values), 
                 pheno_df.identity.values[0], 
                 pheno_df.coverage.values[0], 
                 pheno_df.grade.values[0],
                 pheno_df.databases.values[0]],
                 
            )

            return output

        resfinder_df[["contigs", "identity", "coverage", "grade", "databases"]] = resfinder_df.antibiotic.apply(
            contig_id_info
        )

        return resfinder_df

    def isolate_summary(self):
        """
        Provide an overview of the isolate's results.
        """
        if self.has_json_data is False:
            return self.empty_isolate_summary

        summary_df = pd.DataFrame(
            [
                [
                    self.passport.databases,
                    self.passport.provided_species,
                    self.passport.result_summary,
                ]
            ],
            columns=["databases", "provided_species", "result_summary"],
        )

        summary_df = self.add_analysis_columns(summary_df)

        return summary_df

    @property
    def empty_isolate_summary(self):
        return pd.DataFrame(
            [
                [
                    self.isolate_id,
                    self.time_now,
                    self.passport.key,
                    self.passport.databases,
                    "",
                    "No results found.",
                ]
            ],
            columns=[
                "isolate_id",
                "analysis_date",
                "version",
                "databases",
                "provided_species",
                "result_summary",
            ],
        )


class ResfinderCollector:
    def __init__(self, RESFINDER_dir: str, output_dir: Optional[str] = None,
                 exclude_databases: Optional[List[str]] = None,
    ):
        self.RESFINDER_dir = RESFINDER_dir
        self.output_dir = output_dir
        self._parsers: List[ResfinderParser] = []
        self.isolate_dirs = self.resfinder_dirs()
        self.databases_to_exclude = []
        if exclude_databases is not None:
            self.databases_to_exclude = [db.lower() for db in exclude_databases]

        self.logger = logging.getLogger(__name__)
        # set level to info to get the info messages
        self.logger.setLevel(logging.INFO)
        # log to console
        self.logger.addHandler(logging.StreamHandler())
        self.logger.info(
            f"Found {len(self.isolate_dirs)} isolate directories in {self.RESFINDER_dir}"
        )
        self.logger.info(f"Excluded databases: {', '.join(self.databases_to_exclude) if self.databases_to_exclude else 'None'}")

    def __iter__(self):
        return iter(self._parsers)

    def __len__(self):
        return len(self._parsers)

    def _get_parser(self, isolate_id: str) -> Optional[ResfinderParser]:
        for parser in self._parsers:
            if parser.isolate_id == isolate_id:
                return parser
        return None

    def _add_parser(self, parser: ResfinderParser):
        self._parsers.append(parser)

    def _is_isolate_dir(self, dirpath: str):
        if dirpath.startswith("."):
            return False

        if os.path.isdir(os.path.join(self.RESFINDER_dir, dirpath)) is False:
            return False

        return True

    def resfinder_dirs(self):
        isolate_dirs = os.listdir(self.RESFINDER_dir)
        isolate_dirs = [x for x in isolate_dirs if self._is_isolate_dir(x)]
        return isolate_dirs

    def genes_affected(self, isolate_phenotypes: List[IsolatePhenotypes]):
        all_genes = [
            gene for isolate in isolate_phenotypes for gene in isolate.all_genes()
        ]
        all_genes = list(set(all_genes))

        data_to_keep = []

        for isolate in isolate_phenotypes:
            isolate_line = [isolate.isolate_id, isolate.result_summary]
            isolate.collect_all_genes_affected()
            for gene in all_genes:
                isolate_line.append(isolate.gene_affected(gene))

            data_to_keep.append(isolate_line)

        df = pd.DataFrame(
            data_to_keep, columns=["isolate_id", "result_summary"] + all_genes
        )

        return df

    def variation_summary(
        self, variation_results: pd.DataFrame
    ) -> pd.DataFrame:

        if variation_results.empty:
            return pd.DataFrame(columns=["isolate_id", "Mutation", "Resistance"])

        isolate_ids = variation_results["isolate_id"].unique()

        pointfinder_known = variation_results[
            variation_results["Resistance"] != "Unknown"
        ]

        groups = []
        # for _, group in pointfinder_results.groupby("isolate_id"):
        for isolate_id in isolate_ids:
            group = pointfinder_known[pointfinder_known["isolate_id"] == isolate_id]

            ## merge rows with the same mutation, but different resistance. concatenate the resistance
            group = group.groupby(["isolate_id", "Mutation"]).agg(
                {"Resistance": lambda x: "; ".join(x)}
            )
            group = group.reset_index()

            matrix_df = group.pivot(
                index="isolate_id", columns="Mutation", values="Resistance"
            )

            matrix_df = matrix_df.fillna("")
            matrix_df = matrix_df.reset_index()

            groups.append(matrix_df)

        matrix_df = pd.concat(groups, axis=0)

        return matrix_df
    
    def inspect_directories(self):
        for isolate_dir in self.isolate_dirs:
            if isolate_dir.startswith("."):
                continue

            if os.path.isdir(os.path.join(self.RESFINDER_dir, isolate_dir)) is False:
                continue

            isolate_parser = ResfinderParser(self.RESFINDER_dir, isolate_dir, exclude_databases=self.databases_to_exclude)

            self._add_parser(isolate_parser)

    def collect_all_results(self):

        isolate_summaries = []
        variation_results = []
        region_results = []
        isolate_phenotypes = []

        for isolate_parser in self._parsers:

            isolate_summaries.append(isolate_parser.isolate_summary())

            if isolate_parser.has_json_data is False:
                self.logger.info(
                    f"No results found for isolate {isolate_parser.isolate_id}"
                )
                continue

            isolate_phenotype_results = isolate_parser.collect_phenotype_results()

            isolate_phenotype_results = isolate_parser.extend_seq_region_results(
                isolate_phenotype_results
            )

            isolate_variation_results = isolate_parser.collect_variations_results()

            variation_results.append(isolate_variation_results)
            region_results.append(isolate_phenotype_results)
            isolate_phenotypes.append(isolate_parser.phenotypes)


        region_results = pd.concat(region_results, axis=0)
        isolate_summaries = pd.concat(isolate_summaries, axis=0)
        genes_affected = self.genes_affected(isolate_phenotypes)
        combined_presence_absence = genes_affected.copy()

        if len(variation_results) > 0:
            variation_results = pd.concat(variation_results, axis=0)

            variation_results_summary = self.variation_summary(variation_results)

            combined_presence_absence = pd.merge(
                genes_affected,
                variation_results_summary,
                on="isolate_id",
                how="outer",
            )
        else:
            variation_results = pd.DataFrame()

        return (
            region_results,
            variation_results,
            isolate_summaries,
            combined_presence_absence,
        )

    def collect(self):
        self.logger.info(f"Collecting results from {len(self.isolate_dirs)} isolates..")

        if self.output_dir is None:
            self.output_dir = self.RESFINDER_dir
        else:
            if os.path.isdir(self.output_dir) is False:
                os.mkdir(self.output_dir)

        self.inspect_directories()

        (
            region_results,
            variation_results,
            isolate_summaries,
            combined_presence_absence,
        ) = self.collect_all_results()

        if region_results is None:
            self.logger.info("No results found.")
            return

        region_results.to_csv(
            os.path.join(self.output_dir, "region_results.tsv"),
            sep="\t",
            index=False,
        )
        variation_results.to_csv(
            os.path.join(self.output_dir, "variation_results.tsv"),
            sep="\t",
            index=False,
        )
        isolate_summaries.to_csv(   
            os.path.join(self.output_dir, "isolate_summaries.tsv"),
            sep="\t",
            index=False,
        )

        combined_presence_absence.to_csv(
            os.path.join(self.output_dir, "combined_presence_absence.tsv"),
            sep="\t",
            index=False,
        )

        self.logger.info(f"Results written to {self.output_dir}")
