from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class AssemblyResult:
    contigs: Path | None
    logs: str


@dataclass
class AssemblyStats:
    n_contigs: int
    total_len: int
    n50: int
    largest: int
    gc_pct: float
    contigs_500: int
    completeness: float
    qc_status: str
    qc_flags: list[str]
    misassemblies: int | None = None
    misassembled_contigs: int | None = None
    nga50: int | None = None
    duplication_ratio: float | None = None
    mismatches_per_100kbp: float | None = None
    indels_per_100kbp: float | None = None
    core_genes_found: int | None = None
    core_genes_total: int | None = None

    def to_dict(self) -> dict:
        return {
            "n_contigs":             self.n_contigs,
            "total_len":             self.total_len,
            "n50":                   self.n50,
            "largest":               self.largest,
            "gc_pct":                self.gc_pct,
            "contigs_500":           self.contigs_500,
            "completeness":          self.completeness,
            "misassemblies":         self.misassemblies,
            "misassembled_contigs":  self.misassembled_contigs,
            "nga50":                 self.nga50,
            "duplication_ratio":     self.duplication_ratio,
            "mismatches_per_100kbp": self.mismatches_per_100kbp,
            "indels_per_100kbp":     self.indels_per_100kbp,
            "core_genes_found":      self.core_genes_found,
            "core_genes_total":      self.core_genes_total,
            "qc": {"status": self.qc_status, "flags": self.qc_flags},
        }

    def to_dataframe(self):
        import pandas as pd
        rows = [
            ("n_contigs",    self.n_contigs,    True),
            ("total_len",    self.total_len,    True),
            ("n50",          self.n50,          True),
            ("gc_pct",       self.gc_pct,       True),
            ("completeness", self.completeness, True),
            ("nga50",        self.nga50,        self.nga50 is not None),
        ]
        return pd.DataFrame(
            [(k, v) for k, v, present in rows if present],
            columns=["metric", "value"],
        )


@dataclass
class AMRResult:
    sample_id: str
    chromosomal: dict
    synonymous_variants: dict
    plasmid: dict
    cdc_phenotypes: list
    resistance_category: str
    n_resistance_classes: int
    failure_probability: float
    therapy: dict
    who_matches: list
    mlst: dict
    mosaic_pena: dict
    ngstar: dict
    essential_gene_mutations: dict = None
    essential_gene_synonymous: dict = None

    def to_dict(self) -> dict:
        return {
            "chromosomal":               self.chromosomal,
            "synonymous_variants":       self.synonymous_variants,
            "plasmid":                   self.plasmid,
            "cdc_phenotypes":            self.cdc_phenotypes,
            "resistance_category":       self.resistance_category,
            "n_resistance_classes":      self.n_resistance_classes,
            "failure_probability":       self.failure_probability,
            "therapy":                   self.therapy,
            "who_matches":               self.who_matches,
            "mlst":                      self.mlst,
            "mosaic_pena":               self.mosaic_pena,
            "ngstar":                    self.ngstar,
            "essential_gene_mutations":  self.essential_gene_mutations,
            "essential_gene_synonymous": self.essential_gene_synonymous,
        }

    def to_dataframe(self):
        import pandas as pd
        rows = []
        for gene, muts in self.chromosomal.items():
            for mut in muts:
                rows.append({"gene": gene, "mutation": mut, "source": "chromosomal"})
        for gene, status in self.plasmid.items():
            if status == "present":
                rows.append({"gene": gene, "mutation": None, "source": "plasmid"})
        return pd.DataFrame(rows, columns=["gene", "mutation", "source"])

    @classmethod
    def from_dict(cls, d: dict, sample_id: str = "") -> AMRResult:
        return cls(
            sample_id=sample_id,
            chromosomal=d.get("chromosomal", {}),
            synonymous_variants=d.get("synonymous_variants", {}),
            plasmid=d.get("plasmid", {}),
            cdc_phenotypes=d.get("cdc_phenotypes", []),
            resistance_category=d.get("resistance_category", ""),
            n_resistance_classes=d.get("n_resistance_classes", 0),
            failure_probability=d.get("failure_probability", 0.0),
            therapy=d.get("therapy", {}),
            who_matches=d.get("who_matches", []),
            mlst=d.get("mlst", {}),
            mosaic_pena=d.get("mosaic_pena", {}),
            ngstar=d.get("ngstar", {}),
            essential_gene_mutations=d.get("essential_gene_mutations", {}),
            essential_gene_synonymous=d.get("essential_gene_synonymous", {}),
        )
