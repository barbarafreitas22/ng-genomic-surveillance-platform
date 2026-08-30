_NG_PCT_CRITICAL = 50.0


def detect_alerts(samples: dict) -> list[dict]:
    from .models import AMRResult as _AMRResult
    alerts: list[dict] = []

    for sid, data in samples.items():
        mr = data.get("amr")
        qr = data.get("qc") or {}
        ar = data.get("assembly") or {}

        if not qr.get("error"):
            k2 = qr.get("kraken2") or {}
            ng_pct         = k2.get("ng_pct")
            species_status = k2.get("species_status")

            if species_status == "contaminated":
                ng_str = f"{ng_pct:.1f}%" if ng_pct is not None else "?"
                sev    = "critical" if (ng_pct is not None and ng_pct < _NG_PCT_CRITICAL) else "error"
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "species_contamination",
                    "severity":   sev,
                    "message":    f"Contamination: NG {ng_str}: {sid}",
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
                    "message":    f"Assembly QC fail: {detail}: {sid}",
                })

        if isinstance(mr, _AMRResult):
            res_cat = mr.resistance_category

            if res_cat == "XDR":
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "xdr",
                    "severity":   "critical",
                    "message":    f"XDR profile: {sid}",
                })
            elif res_cat == "MDR":
                alerts.append({
                    "sample_id":  sid,
                    "alert_type": "mdr",
                    "severity":   "critical",
                    "message":    f"MDR profile: {sid}",
                })

    return alerts
