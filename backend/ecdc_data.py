SOURCE = {
    "title": "Confirmed gonorrhoea cases and rates per 100 000 population by country and year, EU/EEA, 2020–2024",
    "publisher": "European Centre for Disease Prevention and Control (ECDC)",
    "report": "Annual Epidemiological Report for 2024",
    "published": "Stockholm: ECDC; 2026",
    "url": "https://www.ecdc.europa.eu/sites/default/files/documents/AER%20gonorrhea%202024.pdf",
    "retrieved": "2026-06-29",
    "note": "NDR = No Data Reported. NRC = No Rate Calculated (sentinel systems or no reporting).",
}

# Notes:
#   NDR = No Data Reported
#   NRC = No Rate Calculated (sentinel surveillance system
#         Belgium, France, Netherlands report sentinel data only;
#         Austria and Germany did not report during this period)
#   ISO 3166-1 numeric codes used to key D3/TopoJSON choropleth.

GONORRHOEA_DATA = {
    40: {
        "name": "Austria",
        "iso2": "AT",
        "no_data": True,
        "sentinel": False,
        "years": {
            2020: {"cases": None, "rate": None},
            2021: {"cases": None, "rate": None},
            2022: {"cases": None, "rate": None},
            2023: {"cases": None, "rate": None},
            2024: {"cases": None, "rate": None},
        },
    },
    56: {
        "name": "Belgium",
        "iso2": "BE",
        "no_data": False,
        "sentinel": True,
        "years": {
            2020: {"cases": 1707,  "rate": None},
            2021: {"cases": 3964,  "rate": None},
            2022: {"cases": 4523,  "rate": None},
            2023: {"cases": 6622,  "rate": None},
            2024: {"cases": 7623,  "rate": None},
        },
    },
    100: {
        "name": "Bulgaria",
        "iso2": "BG",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 17,  "rate": 0.3},
            2021: {"cases": 3,   "rate": 0.0},
            2022: {"cases": 23,  "rate": 0.4},
            2023: {"cases": 56,  "rate": 0.9},
            2024: {"cases": 107, "rate": 1.7},
        },
    },
    191: {
        "name": "Croatia",
        "iso2": "HR",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 13, "rate": 0.3},
            2021: {"cases": 17, "rate": 0.4},
            2022: {"cases": 21, "rate": 0.5},
            2023: {"cases": 26, "rate": 0.7},
            2024: {"cases": 33, "rate": 0.9},
        },
    },
    196: {
        "name": "Cyprus",
        "iso2": "CY",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 7,  "rate": 0.8},
            2021: {"cases": 5,  "rate": 0.6},
            2022: {"cases": 13, "rate": 1.4},
            2023: {"cases": 25, "rate": 2.6},
            2024: {"cases": 21, "rate": 2.2},
        },
    },
    203: {
        "name": "Czechia",
        "iso2": "CZ",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 1672, "rate": 15.6},
            2021: {"cases": 1829, "rate": 17.4},
            2022: {"cases": 2058, "rate": 19.6},
            2023: {"cases": 2593, "rate": 23.9},
            2024: {"cases": 2382, "rate": 21.9},
        },
    },
    208: {
        "name": "Denmark",
        "iso2": "DK",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 3550, "rate": 61.0},
            2021: {"cases": 3611, "rate": 61.8},
            2022: {"cases": 5259, "rate": 89.5},
            2023: {"cases": 5579, "rate": 94.0},
            2024: {"cases": 5091, "rate": 85.4},
        },
    },
    233: {
        "name": "Estonia",
        "iso2": "EE",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 22,  "rate": 1.7},
            2021: {"cases": 54,  "rate": 4.1},
            2022: {"cases": 117, "rate": 8.8},
            2023: {"cases": 137, "rate": 10.0},
            2024: {"cases": 111, "rate": 8.1},
        },
    },
    246: {
        "name": "Finland",
        "iso2": "FI",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 482,  "rate": 8.7},
            2021: {"cases": 510,  "rate": 9.2},
            2022: {"cases": 960,  "rate": 17.3},
            2023: {"cases": 1329, "rate": 23.9},
            2024: {"cases": 1859, "rate": 33.2},
        },
    },
    250: {
        "name": "France",
        "iso2": "FR",
        "no_data": False,
        "sentinel": True,
        "years": {
            2020: {"cases": 5398,  "rate": None},
            2021: {"cases": 7077,  "rate": None},
            2022: {"cases": 8704,  "rate": None},
            2023: {"cases": 10723, "rate": None},
            2024: {"cases": 13533, "rate": None},
        },
    },
    276: {
        "name": "Germany",
        "iso2": "DE",
        "no_data": True,
        "sentinel": False,
        "years": {
            2020: {"cases": None, "rate": None},
            2021: {"cases": None, "rate": None},
            2022: {"cases": None, "rate": None},
            2023: {"cases": None, "rate": None},
            2024: {"cases": None, "rate": None},
        },
    },
    300: {
        "name": "Greece",
        "iso2": "GR",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 161, "rate": 1.5},
            2021: {"cases": 246, "rate": 2.3},
            2022: {"cases": 360, "rate": 3.4},
            2023: {"cases": 457, "rate": 4.4},
            2024: {"cases": 420, "rate": 4.0},
        },
    },
    348: {
        "name": "Hungary",
        "iso2": "HU",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 1261, "rate": 13.0},
            2021: {"cases": 1309, "rate": 13.6},
            2022: {"cases": 1156, "rate": 12.0},
            2023: {"cases": 1345, "rate": 14.0},
            2024: {"cases": 1471, "rate": 15.3},
        },
    },
    352: {
        "name": "Iceland",
        "iso2": "IS",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 93,  "rate": 25.5},
            2021: {"cases": 105, "rate": 28.5},
            2022: {"cases": 158, "rate": 42.0},
            2023: {"cases": 333, "rate": 85.9},
            2024: {"cases": 338, "rate": 88.1},
        },
    },
    372: {
        "name": "Ireland",
        "iso2": "IE",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 2061, "rate": 41.1},
            2021: {"cases": 2349, "rate": 46.4},
            2022: {"cases": 4172, "rate": 80.9},
            2023: {"cases": 6598, "rate": 125.2},
            2024: {"cases": 5832, "rate": 109.0},
        },
    },
    380: {
        "name": "Italy",
        "iso2": "IT",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 333,  "rate": 0.6},
            2021: {"cases": 849,  "rate": 1.4},
            2022: {"cases": 1953, "rate": 3.3},
            2023: {"cases": 2355, "rate": 4.0},
            2024: {"cases": 3105, "rate": 5.3},
        },
    },
    428: {
        "name": "Latvia",
        "iso2": "LV",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 109, "rate": 5.7},
            2021: {"cases": 70,  "rate": 3.7},
            2022: {"cases": 158, "rate": 8.4},
            2023: {"cases": 145, "rate": 7.7},
            2024: {"cases": 157, "rate": 8.4},
        },
    },
    438: {
        "name": "Liechtenstein",
        "iso2": "LI",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 4,  "rate": 10.3},
            2021: {"cases": 5,  "rate": 12.8},
            2022: {"cases": 10, "rate": 25.4},
            2023: {"cases": 10, "rate": 25.2},
            2024: {"cases": 6,  "rate": 15.0},
        },
    },
    440: {
        "name": "Lithuania",
        "iso2": "LT",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 31, "rate": 1.1},
            2021: {"cases": 30, "rate": 1.1},
            2022: {"cases": 38, "rate": 1.4},
            2023: {"cases": 38, "rate": 1.3},
            2024: {"cases": 77, "rate": 2.7},
        },
    },
    442: {
        "name": "Luxembourg",
        "iso2": "LU",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 311, "rate": 49.7},
            2021: {"cases": 417, "rate": 65.7},
            2022: {"cases": 475, "rate": 73.6},
            2023: {"cases": 606, "rate": 91.7},
            2024: {"cases": 582, "rate": 86.6},
        },
    },
    470: {
        "name": "Malta",
        "iso2": "MT",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 94,  "rate": 18.3},
            2021: {"cases": 240, "rate": 46.5},
            2022: {"cases": 228, "rate": 43.8},
            2023: {"cases": 407, "rate": 75.1},
            2024: {"cases": 505, "rate": 89.6},
        },
    },
    528: {
        "name": "Netherlands",
        "iso2": "NL",
        "no_data": False,
        "sentinel": True,
        "years": {
            2020: {"cases": 6826,  "rate": None},
            2021: {"cases": 7966,  "rate": None},
            2022: {"cases": 10601, "rate": None},
            2023: {"cases": 13853, "rate": None},
            2024: {"cases": 13952, "rate": None},
        },
    },
    578: {
        "name": "Norway",
        "iso2": "NO",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 1045, "rate": 19.5},
            2021: {"cases": 555,  "rate": 10.3},
            2022: {"cases": 1858, "rate": 34.2},
            2023: {"cases": 2985, "rate": 54.4},
            2024: {"cases": 3150, "rate": 56.8},
        },
    },
    616: {
        "name": "Poland",
        "iso2": "PL",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 246,  "rate": 0.6},
            2021: {"cases": 287,  "rate": 0.8},
            2022: {"cases": 556,  "rate": 1.5},
            2023: {"cases": 1209, "rate": 3.3},
            2024: {"cases": 1051, "rate": 2.9},
        },
    },
    620: {
        "name": "Portugal",
        "iso2": "PT",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 1068, "rate": 10.3},
            2021: {"cases": 1253, "rate": 12.1},
            2022: {"cases": 2402, "rate": 23.0},
            2023: {"cases": 2768, "rate": 26.3},
            2024: {"cases": 2764, "rate": 26.0},
        },
    },
    642: {
        "name": "Romania",
        "iso2": "RO",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 10, "rate": 0.1},
            2021: {"cases": 22, "rate": 0.1},
            2022: {"cases": 23, "rate": 0.1},
            2023: {"cases": 30, "rate": 0.2},
            2024: {"cases": 38, "rate": 0.2},
        },
    },
    703: {
        "name": "Slovakia",
        "iso2": "SK",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 319, "rate": 5.8},
            2021: {"cases": 414, "rate": 7.6},
            2022: {"cases": 394, "rate": 7.2},
            2023: {"cases": 462, "rate": 8.5},
            2024: {"cases": 412, "rate": 7.6},
        },
    },
    705: {
        "name": "Slovenia",
        "iso2": "SI",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 213, "rate": 10.2},
            2021: {"cases": 292, "rate": 13.8},
            2022: {"cases": 333, "rate": 15.8},
            2023: {"cases": 276, "rate": 13.0},
            2024: {"cases": 194, "rate": 9.1},
        },
    },
    724: {
        "name": "Spain",
        "iso2": "ES",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 10306, "rate": 21.8},
            2021: {"cases": 14605, "rate": 30.8},
            2022: {"cases": 25157, "rate": 53.0},
            2023: {"cases": 34071, "rate": 70.9},
            2024: {"cases": 37169, "rate": 76.4},
        },
    },
    752: {
        "name": "Sweden",
        "iso2": "SE",
        "no_data": False,
        "sentinel": False,
        "years": {
            2020: {"cases": 2692, "rate": 26.1},
            2021: {"cases": 2693, "rate": 25.9},
            2022: {"cases": 3355, "rate": 32.1},
            2023: {"cases": 4232, "rate": 40.2},
            2024: {"cases": 4348, "rate": 41.2},
        },
    },
}

EU_EEA_TOTALS = {
    2020: {"cases": 40051,  "rate": 9.9},
    2021: {"cases": 50777,  "rate": 12.1},
    2022: {"cases": 75065,  "rate": 19.5},
    2023: {"cases": 99270,  "rate": 25.8},
    2024: {"cases": 106331, "rate": 26.9},
}

YEARS = [2020, 2021, 2022, 2023, 2024]

def get_rate(iso_numeric: int, year: int) -> float | None:
    """Return the notification rate for a country and year, or None."""
    entry = GONORRHOEA_DATA.get(iso_numeric)
    if not entry:
        return None
    return entry["years"].get(year, {}).get("rate")


def get_cases(iso_numeric: int, year: int) -> int | None:
    """Return the number of cases for a country and year, or None."""
    entry = GONORRHOEA_DATA.get(iso_numeric)
    if not entry:
        return None
    return entry["years"].get(year, {}).get("cases")


def get_map_data(year: int) -> list[dict]:
    """
    Return a list of dicts ready for D3 choropleth rendering.
    """
    result = []
    for iso, entry in GONORRHOEA_DATA.items():
        yr = entry["years"].get(year, {})
        result.append({
            "id":       iso,
            "name":     entry["name"],
            "iso2":     entry["iso2"],
            "rate":     yr.get("rate"),
            "cases":    yr.get("cases"),
            "sentinel": entry.get("sentinel", False),
            "no_data":  entry.get("no_data", False),
        })
    return result


def get_trend(iso_numeric: int) -> list[dict]:
    """Return year-by-year trend for a country (for sparklines/charts)."""
    entry = GONORRHOEA_DATA.get(iso_numeric)
    if not entry:
        return []
    return [
        {"year": y, "rate": entry["years"][y]["rate"], "cases": entry["years"][y]["cases"]}
        for y in YEARS
    ]


if __name__ == "__main__":
    import json
    print(" ECDC Gonorrhoea Data, EU/EEA 2020-2024")
    print(f"Source: {SOURCE['url']}\n")
    print(f"{'Country':<20} {'2020':>8} {'2021':>8} {'2022':>8} {'2023':>8} {'2024':>8}")
    print("-" * 60)
    for iso, entry in sorted(GONORRHOEA_DATA.items(), key=lambda x: x[1]["name"]):
        rates = []
        for y in YEARS:
            r = entry["years"][y]["rate"]
            rates.append(f"{r:>8.1f}" if r is not None else "     NRC")
        print(f"{entry['name']:<20} {''.join(rates)}")
    print("-" * 60)
    eu = EU_EEA_TOTALS
    print(f"{'EU/EEA total':<20} {eu[2020]['rate']:>8.1f} {eu[2021]['rate']:>8.1f} {eu[2022]['rate']:>8.1f} {eu[2023]['rate']:>8.1f} {eu[2024]['rate']:>8.1f}")
    print(f"\nMap data for 2024 (first 3 entries):")
    print(json.dumps(get_map_data(2024)[:3], indent=2))