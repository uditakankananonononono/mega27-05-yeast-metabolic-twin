import pandas as pd

sw = pd.read_csv('data/raw/swissprot.tsv', sep='\t'); a2n = {}
for _, x in sw.iterrows():
    toks = str(x['gene_id']).split()
    for a in toks: a2n.setdefault(a, toks[0])

s = pd.read_csv('results/strain_design_obligatory_singles.csv')
pos = s[s.obligatory_succinate > 1e-6].copy()
pos['name'] = [a2n.get(g, g) for g in pos.gene]
with open('paper/tab_strain_singles.tex', 'w') as f:
    f.write("\\begin{longtable}{llrr}\n\\caption{Single deletions with obligatory succinate production (min flux at max growth, minimal glucose medium). All other 1{,}127 viable single deletions score zero.}\\label{tab:singles}\\\\\n\\toprule\nORF & gene & growth & obligatory succinate\\\\\n\\midrule\n")
    for _, r in pos.iterrows():
        f.write(f"{r.gene} & {r['name']} & {r.growth:.4f} & {r.obligatory_succinate:.4f}\\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

f5 = pd.read_csv('results/strain_design_5th_ko.csv').head(15).copy()
f5['name'] = [a2n.get(g, g) for g in f5.gene]
with open('paper/tab_strain_5th.tex', 'w') as f:
    f.write("\\begin{longtable}{llrrr}\n\\caption{Top 5th deletions on the Raab quadruple background (baseline obligatory flux 0.0598, growth 0.779).}\\label{tab:fifth}\\\\\n\\toprule\nORF & gene & growth & obligatory & gain\\\\\n\\midrule\n")
    for _, r in f5.iterrows():
        f.write(f"{r.gene} & {r['name']} & {r.growth:.4f} & {r.obligatory_succinate:.4f} & {r.gain_over_quadruple:+.4f}\\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

p = pd.read_csv('results/strain_design_pairs.csv')
bp = p[p.obligatory_succinate > 0.0598 + 1e-9].copy()
bp['n1'] = [a2n.get(g, g) for g in bp.ko1]; bp['n2'] = [a2n.get(g, g) for g in bp.ko2]
with open('paper/tab_strain_pairs.tex', 'w') as f:
    f.write("\\begin{longtable}{llrr}\n\\caption{Double deletions among 23 TCA/glyoxylate target genes exceeding the Raab quadruple baseline (0.0598).}\\label{tab:pairs}\\\\\n\\toprule\nKO1 & KO2 & growth & obligatory succinate\\\\\n\\midrule\n")
    for _, r in bp.iterrows():
        f.write(f"{r.n1} & {r.n2} & {r.growth:.4f} & {r.obligatory_succinate:.4f}\\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")
print('strain tables:', len(pos), len(f5), len(bp))
