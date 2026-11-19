"""Command-line interface for the complete Python pipeline."""

import argparse
import json
import sys
from pathlib import Path

from .digest import FORWARD_27F, REVERSE_1492R, generate_library
from .gel import LADDERS, analyze_gel
from .lanes import read_lane_map
from .match import (
    identify_samples,
    identify_samples_consensus,
    validate_controls,
    write_dendrograms,
)
from .profiles import read_profiles
from .report import write_html_report


def _list(value, cast=str):
    if value is None or isinstance(value, (list, tuple)):
        return value
    return [cast(item.strip()) for item in str(value).split(",") if item.strip()]


def _crop(value):
    values = _list(value, int)
    if value is not None and len(values) != 4:
        raise argparse.ArgumentTypeError("crop must be x,y,width,height")
    return tuple(values) if values is not None else None


def _add_gel_options(parser, configurable=False):
    required = not configurable
    parser.add_argument("--image", required=required, help="Gel image (PNG/JPEG/TIFF)")
    parser.add_argument("--enzyme", required=required, help="Restriction enzyme used for the gel")
    parser.add_argument("--n-lanes", type=int, help="Total lanes, including ladder; inferred from --lane-map")
    parser.add_argument("--lane-map", help="CSV with lane,role,sample,sample_type,expected_match")
    parser.add_argument("--output-dir", default="results/gel" if not configurable else None)
    parser.add_argument("--ladder", choices=sorted(LADDERS), help="Built-in ladder preset")
    parser.add_argument("--ladder-sizes", help="Comma-separated bp sizes, largest to smallest")
    parser.add_argument("--ladder-lane", type=int, help="One-based primary ladder lane number")
    parser.add_argument("--secondary-ladder-lane", type=int, help="Optional opposite ladder for smile correction")
    parser.add_argument("--ladder-y", help="Optional comma-separated primary ladder pixel y positions")
    parser.add_argument("--secondary-ladder-y", help="Optional comma-separated secondary ladder pixel y positions")
    parser.add_argument("--lane-centers", help="Optional comma-separated lane x centers")
    parser.add_argument("--sample-names", help="Sample names, excluding the ladder")
    parser.add_argument("--crop", type=_crop, help="Gel crop as x,y,width,height")
    parser.add_argument("--polarity", choices=["bright", "dark"], default="bright")
    parser.add_argument("--prominence", type=float, default=2.5, help="Band detection sensitivity; lower detects more")
    parser.add_argument("--min-band-distance", type=int, help="Minimum peak spacing in pixels; controls close doublets")
    parser.add_argument("--experiment-metadata", help="JSON file for gel percentage, voltage, run time, operator, notes")


def _add_digest_options(parser, configurable=False):
    parser.add_argument("--fasta", required=not configurable, help="Reference 16S FASTA")
    parser.add_argument("--enzymes", default="HaeIII,HhaI")
    parser.add_argument("--forward-primer", default=FORWARD_27F)
    parser.add_argument("--reverse-primer", default=REVERSE_1492R)
    parser.add_argument("--primer-mismatches", type=int, default=2)
    parser.add_argument("--allow-unclipped", action="store_true", help="Digest available 16S sequence when primer sites are absent; terminal fragments are approximate")
    parser.add_argument("--double-digest", action="append", default=[], help="Enzyme pair, e.g. HaeIII+HhaI; repeatable")


def _add_match_options(parser):
    parser.add_argument("--tolerance", type=float, default=20, help="Band matching tolerance in bp")
    parser.add_argument("--min-bp", type=float, default=100)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--exclude-sample", action="append", default=[], help="Sample/marker to exclude; repeatable")
    parser.add_argument("--combine-enzymes", action="store_true", help="Consensus ranking across all sample enzymes")


def _load_run_defaults(argv):
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config")
    known, _ = pre.parse_known_args(argv)
    if not known.config:
        return {}
    config = json.loads(Path(known.config).read_text())
    defaults = {key: value for key, value in config.items() if not isinstance(value, dict)}
    for section in ("gel", "digest", "matching"):
        defaults.update(config.get(section, {}))
    return defaults


def build_parser(defaults=None):
    parser = argparse.ArgumentParser(prog="ardronics", description="Gel-to-16S ARDRA analysis")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gel = subparsers.add_parser("gel", help="Read bands and molecular sizes from a gel")
    _add_gel_options(gel)

    library = subparsers.add_parser("library", help="Generate an in-silico digest library")
    _add_digest_options(library)
    library.add_argument("--output", default="results/library.csv")

    identify = subparsers.add_parser("identify", help="Match sample profiles to a library")
    identify.add_argument("--samples", required=True)
    identify.add_argument("--library", required=True)
    identify.add_argument("--output-dir", default="results/identification")
    _add_match_options(identify)

    run = subparsers.add_parser("run", help="Run gel reading, library generation, and identification")
    run.add_argument("--config", help="JSON configuration; explicit CLI options take precedence")
    _add_gel_options(run, configurable=True)
    _add_digest_options(run, configurable=True)
    _add_match_options(run)
    run.set_defaults(**(defaults or {}))
    return parser


def _gel_arguments(args):
    sizes = _list(args.ladder_sizes, int)
    if sizes is None and args.ladder:
        sizes = LADDERS[args.ladder]
    arguments = {
        "image_path": args.image, "enzyme": args.enzyme, "output_dir": args.output_dir,
        "n_lanes": args.n_lanes, "ladder_lane": args.ladder_lane,
        "secondary_ladder_lane": args.secondary_ladder_lane, "ladder_sizes": sizes,
        "ladder_y": _list(args.ladder_y, float), "secondary_ladder_y": _list(args.secondary_ladder_y, float),
        "sample_names": _list(args.sample_names),
        "crop": _crop(args.crop), "lane_centers": _list(args.lane_centers, int),
        "polarity": args.polarity, "prominence": args.prominence,
        "min_band_distance": args.min_band_distance,
        "experiment_metadata": json.loads(Path(args.experiment_metadata).read_text()) if args.experiment_metadata else None,
    }
    args.control_expectations = {}
    if args.lane_map:
        if args.n_lanes is not None or args.sample_names is not None or args.ladder_lane is not None or args.secondary_ladder_lane is not None:
            raise ValueError("Lane map is authoritative: do not also supply n-lanes, sample-names, or ladder-lane options")
        lane_map = read_lane_map(args.lane_map)
        arguments.update({
            "n_lanes": lane_map.n_lanes, "ladder_lane": lane_map.ladder_lane,
            "secondary_ladder_lane": lane_map.secondary_ladder_lane,
            "sample_names": lane_map.sample_names, "lane_roles": lane_map.roles,
            "sample_types": lane_map.sample_types,
        })
        args.exclude_sample = [*getattr(args, "exclude_sample", []), *lane_map.excluded_samples]
        args.control_expectations = lane_map.expected_matches
    if arguments["n_lanes"] is None:
        raise ValueError("Supply --n-lanes or --lane-map")
    return arguments


def _identify(args, samples_path, library_path, output_dir):
    samples = read_profiles(samples_path)
    samples = samples.loc[~samples["sample"].isin(["ladder", *args.exclude_sample])]
    library = read_profiles(library_path)
    matcher = identify_samples_consensus if args.combine_enzymes else identify_samples
    report = matcher(samples, library, args.tolerance, args.min_bp, args.top_n)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_dir / "identification.csv", index=False, float_format="%.4f")
    write_dendrograms(samples, library, output_dir, args.tolerance, args.min_bp)
    return report


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    try:
        defaults = _load_run_defaults(argv) if argv and argv[0] == "run" else {}
        parser = build_parser(defaults)
        args = parser.parse_args(argv)
        if args.command == "gel":
            result = analyze_gel(**_gel_arguments(args))
            print(f"Detected {len(result.bands)} bands in {len(result.lane_centers)} lanes")
            print(f"Outputs: {result.output_dir}")
            for warning in result.warnings:
                print(f"WARNING: {warning}", file=sys.stderr)
        elif args.command == "library":
            library = generate_library(**_digest_arguments(args, args.output))
            print(f"Wrote {len(library)} fragments for {library['sample'].nunique()} references to {args.output}")
        elif args.command == "identify":
            report = _identify(args, args.samples, args.library, args.output_dir)
            print(f"Wrote {len(report)} ranked matches to {Path(args.output_dir) / 'identification.csv'}")
        elif args.command == "run":
            required = ["image", "enzyme", "fasta"]
            missing = [f"--{name.replace('_', '-')}" for name in required if getattr(args, name, None) is None]
            if args.n_lanes is None and args.lane_map is None:
                missing.append("--n-lanes or --lane-map")
            if missing:
                raise ValueError(f"Missing required run settings: {', '.join(missing)}")
            if not args.output_dir:
                args.output_dir = "results/run"
            if args.ladder_sizes is None and args.ladder is None:
                raise ValueError("run requires --ladder or --ladder-sizes")
            output = Path(args.output_dir)
            gel_args = _gel_arguments(args)
            gel_args["output_dir"] = output / "gel"
            gel = analyze_gel(**gel_args)
            generate_library(**_digest_arguments(args, output / "library.csv"))
            report = _identify(args, output / "gel" / "profiles.csv", output / "library.csv", output)
            if args.control_expectations:
                controls = validate_controls(report, args.control_expectations)
                controls.to_csv(output / "control_validation.csv", index=False)
                if not controls["passed"].all():
                    print("WARNING: One or more positive controls failed; see control_validation.csv", file=sys.stderr)
            report_path = write_html_report(output, report, gel.warnings)
            print(f"Analysis complete: {report_path}")
        return 0
    except (ValueError, FileNotFoundError, json.JSONDecodeError) as error:
        print(f"ardronics: error: {error}", file=sys.stderr)
        return 2


def _digest_arguments(args, output):
    pairs = []
    for pair in args.double_digest:
        enzymes = tuple(pair.split("+"))
        if len(enzymes) != 2 or not all(enzymes):
            raise ValueError(f"Invalid double digest: {pair}")
        pairs.append(enzymes)
    return {
        "fasta_path": args.fasta, "output_csv": output, "enzymes": _list(args.enzymes),
        "forward_primer": args.forward_primer, "reverse_primer": args.reverse_primer,
        "max_primer_mismatches": args.primer_mismatches, "double_digests": pairs,
        "allow_unclipped": args.allow_unclipped,
    }


if __name__ == "__main__":
    raise SystemExit(main())

