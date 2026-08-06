import argparse
import json
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.paths import PROJECT_ROOT

from backend.amr import (
    GENE_DB_CROM,
    _SKIP_CHROM_GENES,
    call_variants,
    parse_vcf,
)


def build_panel(genomes_dir: Path, output_path: Path) -> None:
    genome_files = sorted(
        list(genomes_dir.glob("*.fasta"))
        + list(genomes_dir.glob("*.fa"))
        + list(genomes_dir.glob("*.fna"))
    )

    if not genome_files:
        print(f"No FASTA files found in {genomes_dir}", file=sys.stderr)
        sys.exit(1)

    n = len(genome_files)
    print(f"Reference panel: {n} genome(s) in {genomes_dir}")

    gene_refs = [
        f for f in sorted(GENE_DB_CROM.glob("*.fasta"))
        if f.stem not in _SKIP_CHROM_GENES
    ]
    print(f"Genes scanned: {len(gene_refs)}")

    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for i, genome in enumerate(genome_files):
        print(f"  [{i+1}/{n}] {genome.name}")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            for ref_gene in gene_refs:
                try:
                    vcf = call_variants(ref_gene, genome, tmp_path)
                    _muts, syn_muts = parse_vcf(vcf, ref_gene=ref_gene)
                    for s in syn_muts:
                        counts[ref_gene.stem][s["nuc"]] += 1
                except Exception:
                    pass

    frequencies: dict[str, dict[str, float]] = {}
    for gene, mut_counts in sorted(counts.items()):
        frequencies[gene] = {
            mut: round(count / n, 4)
            for mut, count in sorted(mut_counts.items())
        }

    panel = {
        "meta": {
            "description": "Synonymous mutation frequency panel for N. gonorrhoeae surveillance",
            "n_genomes": n,
            "genomes": [f.name for f in genome_files],
        },
        "frequencies": frequencies,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(panel, fh, indent=2)

    total_muts = sum(len(v) for v in frequencies.values())
    print(f"Done — {total_muts} unique synonymous variants across {n} genome(s)")
    print(f"Panel saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", required=True, type=Path, metavar="DIR",
                        help="Directory containing assembled genome FASTAs (.fasta/.fa/.fna)")
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "data" / "synonymous_freq_panel.json",
                        metavar="FILE", help="Output JSON path (default: data/synonymous_freq_panel.json)")
    args = parser.parse_args()

    if not args.genomes.is_dir():
        print(f"Error: {args.genomes} is not a directory", file=sys.stderr)
        sys.exit(1)

    build_panel(args.genomes, args.output)
