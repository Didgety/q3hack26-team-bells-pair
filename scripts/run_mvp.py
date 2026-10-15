"""Regenerate every MVP table, shot record, figure, and the run manifest."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from q3hack.experiment import Config, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    parser.add_argument("--quick", action="store_true", help="Smaller complete pipeline for notebook/iteration; not the main statistical run")
    args = parser.parse_args()
    config = Config(short_shots=256, calibration_repeats=4, evaluation_shots=512,
                    bootstrap_repeats=200, budgets=(0, 128, 512)) if args.quick else Config()
    run(args.output, config)


if __name__ == "__main__":
    main()
