import os
from io import StringIO

import newick

try:
    from Bio import Phylo as _Phylo
    _HAS_BIOPYTHON = True
except ImportError:
    _HAS_BIOPYTHON = False

OUTGROUPS = ["lactamica_genomic", "mendinitis_genomic"]


def reroot_tree(newick_str: str) -> str:
    if not _HAS_BIOPYTHON or not newick_str:
        return newick_str
    try:
        _PH = "XUSCOREX"
        encoded = newick_str.replace("_", _PH)
        encoded_outgroups = [og.replace("_", _PH) for og in OUTGROUPS]

        tree = _Phylo.read(StringIO(encoded), "newick")
        terminals = {c.name for c in tree.get_terminals() if c.name}
        present = [og for og in encoded_outgroups if og in terminals]

        if len(present) >= 2:
            clades = [tree.find_any(name=og) for og in present]
            mrca = tree.common_ancestor(clades)
            tree.root_with_outgroup(mrca)
        elif len(present) == 1:
            tree.root_with_outgroup({"name": present[0]})
        else:
            tree.root_at_midpoint()

        out = StringIO()
        _Phylo.write(tree, out, "newick")
        return out.getvalue().strip().replace(_PH, "_")
    except Exception:
        return newick_str


def load_tree_newick(tree_path):
    if not tree_path or not os.path.exists(tree_path):
        return None, f"Treefile not found at:\n{tree_path}"
    try:
        with open(tree_path) as f:
            newick_str = f.read().strip()
        if not newick_str:
            return None, "Treefile exists but is empty."
        return newick_str, None
    except Exception as e:
        return None, f"Could not load tree: {e}"


def newick_to_tree_json(newick_str, clusters=None, upload_names=None, backbone_meta=None, resistance=None):
    """
    Convert Newick string to hierarchical JSON for D3 dendrogram.
    Node colours: uploaded = red, outgroup = grey, reference panel = blue.
    backbone_meta: dict mapping node name → metadata dict (from backbone_metadata.json).
    resistance: optional dict mapping node name → "resistant"/"susceptible"/"no_data",
    used to highlight nodes affected by a selected antibiotic.
    """
    clusters      = clusters      or {}
    backbone_meta = backbone_meta or {}
    resistance    = resistance    or {}
    _upload_set   = set(upload_names) if upload_names else set()
    newick_str = reroot_tree(newick_str)
    trees = newick.loads(newick_str)
    root = trees[0]
    counter = [0]

    def _color(node_type):
        if node_type == "internal":
            return "#9ca3af"
        if node_type == "outgroup":
            return "#6b7280"
        if node_type == "leaf_upload":
            return "#dc2626"
        return "#1d4ed8"

    def to_dict(n):
        if n.name:
            node_id = n.name
            if n.name in OUTGROUPS:
                node_type = "outgroup"
            elif n.name in _upload_set:
                node_type = "leaf_upload"
            else:
                node_type = "leaf_base"
        else:
            node_id = f"internal_{counter[0]}"
            counter[0] += 1
            node_type = "internal"

        cluster_id = str(clusters.get(node_id, ""))
        try:
            length = float(n.length) if n.length is not None else 0.0
        except (ValueError, TypeError):
            length = 0.0

        d = {
            "name":       n.name or "",
            "id":         node_id,
            "length":     length,
            "type":       node_type,
            "cluster":    cluster_id,
            "color":      _color(node_type),
            "resistance": resistance.get(node_id),
        }
        if node_id in backbone_meta:
            d["backbone_info"] = backbone_meta[node_id]
        if n.descendants:
            d["children"] = [to_dict(child) for child in n.descendants]
        return d

    return to_dict(root)
