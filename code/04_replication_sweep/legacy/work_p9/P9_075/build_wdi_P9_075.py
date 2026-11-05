# Builder (run once, 2026-10-11) for work_p9/P9_075/wdi_P9_075.csv: World Bank WDI API v2, all countries, 1965:2002.
# Usage: python3 build_wdi_P9_075.py <outdir>   (downloads 4 JSON pages, writes wdi_P9_075.csv: iso3c, year, debt_ppg, gdp, exports, inflation)
import json, csv, sys, urllib.request
IND = ["DT.DOD.DPPG.CD", "NY.GDP.MKTP.CD", "NE.EXP.GNFS.CD", "FP.CPI.TOTL.ZG"]
S = {}
for ind in IND:
    url = f"https://api.worldbank.org/v2/country/all/indicator/{ind}?format=json&per_page=20000&date=1965:2002"
    for r in json.load(urllib.request.urlopen(url))[1]:
        S.setdefault((r["countryiso3code"], r["date"]), {})[ind] = r["value"]
w = csv.writer(open(sys.argv[1] + "/wdi_P9_075.csv", "w", newline=""), lineterminator="\n")
w.writerow(["iso3c", "year", "debt_ppg", "gdp", "exports", "inflation"])
for (iso, y), v in sorted(S.items()):
    if iso: w.writerow([iso, y] + ["" if v.get(i) is None else v.get(i) for i in IND])
