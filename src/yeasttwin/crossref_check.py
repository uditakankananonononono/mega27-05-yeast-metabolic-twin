"""Crossref bibliography-verification lane.

Parses the embedded thebibliography of the paper, queries the Crossref
REST API for each reference, and flags entries whose top hit shows weak
title-token overlap (candidate mis-citations). Paper QA, honestly reported.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
UA = {'User-Agent': 'mega27-yeast-twin/1.0 (mailto:mega27@example.org)'}
STOP = {"the", "of", "and", "a", "in", "for", "with", "on", "by", "to",
        "an", "from", "using", "via", "its"}


def tokens(s):
    return {w for w in re.findall(r"[a-z]{3,}", s.lower()) if w not in STOP}


def crossref_top(query):
    q = urllib.parse.quote(query)
    url = f"https://api.crossref.org/works?query.bibliographic={q}&rows=1"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        items = d["message"]["items"]
        return items[0] if items else {}
    except Exception:
        return {}


def main():
    tex = (ROOT / "paper" / "main.tex").read_text()
    bib = tex.split("\\begin{thebibliography}", 1)[1]
    items = re.findall(r"\\bibitem\{([^}]+)\}\s*(.+?)(?=\\bibitem|\\end\{thebibliography\})",
                       bib, re.S)
    rows = []
    for key, text in items:
        text = re.sub(r"\\emph\{([^}]*)\}", r"\1", text)
        text = re.sub(r"\s+", " ", text).strip()
        top = crossref_top(text)
        title = (top.get("title") or [""])[0]
        ref_tok, hit_tok = tokens(text), tokens(title)
        ov = (len(ref_tok & hit_tok) / len(ref_tok & hit_tok | hit_tok)) if ref_tok else 0.0
        rows.append({"key": key, "reference": text[:160],
                     "crossref_title": title[:160],
                     "doi": top.get("DOI", ""),
                     "token_overlap": round(ov, 3)})
        print(key, round(ov, 2), top.get("DOI", ""), flush=True)
        time.sleep(0.3)
    with open(ROOT / "results" / "crossref_bib_check.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    verified = sum(1 for r in rows if r["token_overlap"] >= 0.2 and r["doi"])
    summary = {"n_references": len(rows), "n_crossref_matched": verified}
    (ROOT / "results" / "crossref_bib_summary.json").write_text(
        json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
