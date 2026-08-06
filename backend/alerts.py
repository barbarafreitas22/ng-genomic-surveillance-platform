_SEQ_COVERAGE_MIN  = 30.0
_SEQ_BREADTH_MIN   = 80.0
_FAILURE_PROB_HIGH = 0.50
_NG_PCT_WARN       = 80.0   
_NG_PCT_CRITICAL   = 50.0  


def detect_alerts(samples: dict) -> list[dict]:
    from .models import AMRResult as _AMRResult
    alerts: list[dict] = []

    for sid, data in samples.items():
        mr = data.get("amr")
        qr = data.get("qc") or {}
        ar = data.get("assembly") or {}

        if isinstance(mr, _AMRResult):
            res_cat = mr.resistance_category
            prob    = mr.failure_probability
            cdc     = mr.cdc_phenotypes or []
            pheno   = "; ".join(
                c.replace("_", " ").capitalize() for c in cdc if c != "wildtype"
            ) or "—"

            if res_cat == "MDR":
                _who = mr.who_matches or []
                _is_xdr = any(m.get("mdr_class") == "XDR" for m in _who)
                if _is_xdr:
                    alerts.append({
                        "sample_id":  sid,
                        "alert_type": "xdr",
                        "severity":   "critical",
                        "message":    f"XDR profile — {sid}",
                    })
                else:
                    alerts.append({
                        "sample_id":  sid,
                        "alert_type": "mdr",
                        "severity":   "critical",
                        "message":    f"MDR profile — {sid}",
                    })
            elif res_cat == "high_resistance":
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "high_resistance",
                    "severity":   "error",
                    "message":    f"High-level resistance — {sid}",
                })

            if prob is not None and prob >= _FAILURE_PROB_HIGH:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "treatment_failure_risk",
                    "severity":   "error",
                    "message":    f"Treatment failure risk {prob*100:.0f}% — {sid}",
                })

            if "ceftriaxone_reduced_susceptibility" in cdc:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "ceftriaxone_resistance",
                    "severity":   "critical",
                    "message":    f"Ceftriaxone reduced susceptibility — {sid}",
                })

            if "high_level_azithromycin_resistance" in cdc:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "hl_azithromycin",
                    "severity":   "error",
                    "message":    f"High-level azithromycin resistance — {sid}",
                })

            if "ciprofloxacin_resistant" in cdc:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "ciprofloxacin_resistant",
                    "severity":   "warning",
                    "message":    f"Ciprofloxacin resistance (gyrA/parC) — {sid}",
                })

            mos = mr.mosaic_pena or {}
            if mos.get("mosaic_suspected"):
                mos_class = mos.get("allele_class") or "unknown"
                mos_id    = mos.get("identity")
                id_str    = f" {mos_id:.1f}%" if mos_id else ""
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "mosaic_pena",
                    "severity":   "error",
                    "message":    f"Mosaic penA {mos_class}{id_str} — {sid}",
                })

            plasmid_res = {
                gene: status
                for gene, status in mr.plasmid.items()
                if status == "present"
            }
            if plasmid_res:
                genes_str = ", ".join(sorted(plasmid_res.keys()))
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "plasmid_resistance",
                    "severity":   "warning",
                    "message":    f"Plasmid genes: {genes_str} — {sid}",
                })

            ngstar = mr.ngstar or {}
            if ngstar.get("novel") and not ngstar.get("incomplete"):
                st = ngstar.get("ST") or "?"
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "novel_ngstar",
                    "severity":   "warning",
                    "message":    f"Novel NG-STAR ST-{st} — {sid}",
                })

        if not ar.get("error"):
            asm_qc = (ar.get("stats") or {}).get("qc") or ar.get("qc") or {}
            if asm_qc.get("status") == "fail":
                flags  = asm_qc.get("flags") or []
                detail = flags[0] if flags else "below threshold"
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "assembly_qc_fail",
                    "severity":   "error",
                    "message":    f"Assembly QC fail: {detail} — {sid}",
                })
            elif asm_qc.get("status") == "warn":
                flags  = asm_qc.get("flags") or []
                detail = flags[0] if flags else ""
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "assembly_qc_warn",
                    "severity":   "warning",
                    "message":    f"Assembly QC warning: {detail} — {sid}" if detail else f"Assembly QC warning — {sid}",
                })

        if not qr.get("error"):
            met  = qr.get("metrics") or {}
            cov  = met.get("mean_coverage")
            b10x = met.get("pct_breadth_10x")
            k2   = qr.get("kraken2") or {}
            ng_pct         = k2.get("ng_pct")
            species_status = k2.get("species_status")

            if cov is not None and cov < _SEQ_COVERAGE_MIN:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "low_coverage",
                    "severity":   "warning",
                    "message":    f"Low coverage {cov:.0f}× — {sid}",
                })

            if b10x is not None and b10x < _SEQ_BREADTH_MIN:
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "low_breadth",
                    "severity":   "warning",
                    "message":    f"Low breadth {b10x:.0f}% at ≥10× — {sid}",
                })

            if species_status == "contaminated":
                ng_str = f"{ng_pct:.1f}%" if ng_pct is not None else "?"
                sev    = "critical" if (ng_pct is not None and ng_pct < _NG_PCT_CRITICAL) else "error"
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "species_contamination",
                    "severity":   sev,
                    "message":    f"Contamination: NG {ng_str} — {sid}",
                })

    return alerts
