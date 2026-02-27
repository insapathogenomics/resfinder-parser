import argparse
from resfinder_parser import ResfinderCollector


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--resfinder_dir",
        help="Path to the directory containing the resfinder results.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="Path to the directory to write the results to.",
        required=False,
    )

    parser.add_argument(
        "-e",
        "--exclude_databases",
        help="Specify databases to exclude. This argument can be used multiple times.",
        action="append",
        required=False,
    )

    args = parser.parse_args()

    excluded_databases = args.exclude_databases if args.exclude_databases else []
    excluded_databases = [db.lower() for db in excluded_databases]
    excluded_databases = list(set(excluded_databases))

    collector = ResfinderCollector(args.resfinder_dir, args.output_dir, exclude_databases=excluded_databases)
    collector.collect()


if __name__ == "__main__":
    main()
