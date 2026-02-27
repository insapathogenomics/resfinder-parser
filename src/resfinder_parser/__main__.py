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
    args = parser.parse_args()

    collector = ResfinderCollector(args.resfinder_dir, args.output_dir)
    collector.collect()


if __name__ == "__main__":
    main()
