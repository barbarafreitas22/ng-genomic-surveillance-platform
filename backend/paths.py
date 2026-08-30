import configparser as _configparser
from pathlib import Path

PROJECT_ROOT      = Path(__file__).resolve().parents[1]
CONFIG_PATH       = PROJECT_ROOT / "config.ini"

DATA_DIR          = PROJECT_ROOT / "data"
REFERENCE_FA1090  = DATA_DIR / "reference" / "FA1090.fasta"

GENE_DB_CROM      = DATA_DIR / "genes" / "chromosomal"
GENE_DB_PLASM     = DATA_DIR / "genes" / "plasmid"
GENE_DB_ESSENTIAL = DATA_DIR / "genes" / "essential"

CORE_GENOME_FASTA = DATA_DIR / "genes" / "core_genome.fasta"

_config = _configparser.ConfigParser()
_config.read(CONFIG_PATH)
_results_dir_raw = _config.get("PATHS", "RESULTS_DIR", fallback="results")
RESULTS_DIR = (
    PROJECT_ROOT / _results_dir_raw
    if not Path(_results_dir_raw).is_absolute()
    else Path(_results_dir_raw)
)

MLST_PROFILES_CACHE   = RESULTS_DIR / "mlst_profiles_cache.tsv"
NGSTAR_PROFILES_CACHE = RESULTS_DIR / "ngstar_profiles_cache.tsv"
NGMAST_PROFILES_CACHE = RESULTS_DIR / "ngmast_profiles_cache.tsv"

PROJECTS_DIR      = PROJECT_ROOT / "app" / "projects"

PHYLOGENY_DIR     = PROJECT_ROOT / "backend" / "phylogeny"
RUNS_DIR          = PHYLOGENY_DIR / "runs"

PHYLOGENY_DATA_DIR      = DATA_DIR / "phylogeny"
BASE_GENOMES_DIR        = PHYLOGENY_DATA_DIR / "genomes"
WHOF_REF                = BASE_GENOMES_DIR / "WHOF_genomic.fna"
BACKBONE_METADATA_PATH  = PHYLOGENY_DATA_DIR / "backbone_metadata.json"

SKA_CACHE_DIR           = PHYLOGENY_DATA_DIR / "ska_cache"
BACKBONE_SKETCH         = SKA_CACHE_DIR / "backbone.skf"
BACKBONE_HASH_FILE      = SKA_CACHE_DIR / "backbone.md5"
BACKBONE_FILE_LIST      = SKA_CACHE_DIR / "backbone_file_list.tsv"
SAMPLE_SKETCH_CACHE_DIR = SKA_CACHE_DIR / "sample_sketch_cache"
SAMPLE_SKETCH_MANIFEST  = SKA_CACHE_DIR / "sample_sketch_manifest.json"
POPPUNK_DB_DIR          = PHYLOGENY_DATA_DIR / "poppunk_db"

PYNGOST_DB_DIR          = PHYLOGENY_DATA_DIR / "pyngost_db"
CGMLST_SCHEMA_DIR       = PHYLOGENY_DATA_DIR / "cgmlst_schema"

SNP_GENOGROUP_THRESHOLD    = 2000
