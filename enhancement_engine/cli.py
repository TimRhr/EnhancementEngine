# Command-line interface for Enhancement Engine
import argparse
import os
from typing import List, Optional

from .core.engine import EnhancementEngine


def _parse_gene_list(value: str) -> List[str]:
    """Parse gene list from comma-separated string or file path."""
    if os.path.isfile(value):
        with open(value, "r", encoding="utf-8") as fh:
            return [line.strip() for line in fh if line.strip()]
    return [v.strip() for v in value.split(",") if v.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Enhancement Engine CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Single gene analysis
    analyze = subparsers.add_parser("analyze", help="Analyze a single gene")
    analyze.add_argument("gene", help="Gene symbol to analyze")
    analyze.add_argument("-e", "--email", required=True, help="Email for NCBI access")
    analyze.add_argument("-v", "--variant", default="enhancement_variant", help="Variant to analyze")
    analyze.add_argument("-t", "--target-tissue", default="general", help="Target tissue")

    # Batch analysis
    batch = subparsers.add_parser("batch", help="Analyze multiple genes")
    batch.add_argument("genes", help="Comma-separated genes or file path")
    batch.add_argument("-e", "--email", required=True, help="Email for NCBI access")
    batch.add_argument("-v", "--variants", help="Comma-separated variants matching genes")

    # Cache management
    cache = subparsers.add_parser("cache", help="Cache management")
    cache.add_argument("-e", "--email", required=True, help="Email for NCBI access")
    cache_sub = cache.add_subparsers(dest="cache_cmd", required=True)
    cache_sub.add_parser("clear", help="Clear cached data")
    cache_sub.add_parser("stats", help="Show cache statistics")

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "analyze":
        engine = EnhancementEngine(args.email)
        report = engine.analyze_gene(args.gene, args.variant, args.target_tissue)
        print(f"Gene: {report.gene_name}")
        print(f"Feasibility score: {report.feasibility_score:.1f}")
        print(f"Safety score: {report.safety_assessment.overall_score:.1f}")
        gain = report.predicted_effect.enhancement_gain.improvement_factor
        print(f"Enhancement gain: {gain:.1f}x")

    elif args.command == "batch":
        genes = _parse_gene_list(args.genes)
        variants = None
        if args.variants:
            variants = [v.strip() for v in args.variants.split(",") if v.strip()]
        engine = EnhancementEngine(args.email)
        batch_report = engine.batch_analysis(genes, variants)
        print(f"Successfully analyzed {batch_report.successful_analyses}/"
              f"{batch_report.total_genes} genes")
        for gene, rep in batch_report.results.items():
            safety = rep.safety_assessment.overall_score
            feas = rep.feasibility_score
            print(f"{gene}: feasibility {feas:.1f}, safety {safety:.1f}")
        if batch_report.failed_analyses:
            print("Failed analyses:")
            for gene, msg in batch_report.failed_analyses.items():
                print(f" - {gene}: {msg}")

    elif args.command == "cache":
        engine = EnhancementEngine(args.email)
        if args.cache_cmd == "clear":
            engine.clear_cache()
            print("Cache cleared")
        elif args.cache_cmd == "stats":
            stats = engine.get_database_stats()
            for key, value in stats.items():
                print(f"{key}: {value}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
