"""
01_fasta_qc.py
Basic FASTA QC
"""

from pathlib import Path
from statistics import mean, median
from collections import Counter
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt

INPUT_FASTA = "genes_nucleotide.fasta"
VALID_BASES = set("ATGCN")

def gc_content(seq:str)->float:
    seq=seq.upper()
    return round(((seq.count("G")+seq.count("C"))/len(seq))*100,2) if seq else 0

def analyse_fasta(fp):
    rows=[]
    seq_seen={}
    dup=[]
    for rec in tqdm(list(SeqIO.parse(fp,"fasta")),desc="Reading"):
        s=str(rec.seq).upper()
        invalid=",".join(sorted(set(s)-VALID_BASES))
        if s in seq_seen:
            dup.append({"Gene_ID":rec.id,"Duplicate_Of":seq_seen[s]})
        else:
            seq_seen[s]=rec.id
        rows.append({
            "Gene_ID":rec.id,
            "Length":len(s),
            "GC(%)":gc_content(s),
            "N_Count":s.count("N"),
            "Invalid_Characters":invalid
        })
    return pd.DataFrame(rows),pd.DataFrame(dup)

def main():
    fasta=Path(INPUT_FASTA)
    if not fasta.exists():
        raise FileNotFoundError(INPUT_FASTA)
    out=Path("output")
    out.mkdir(exist_ok=True)

    df,dup=analyse_fasta(fasta)
    invalid=df[df["Invalid_Characters"]!=""]

    summary=pd.DataFrame({
        "Metric":[
            "Total genes","Total bases","Minimum length","Maximum length",
            "Mean length","Median length","Overall GC (%)",
            "Total ambiguous bases","Invalid sequences","Duplicate sequences"
        ],
        "Value":[
            len(df),
            int(df["Length"].sum()),
            int(df["Length"].min()),
            int(df["Length"].max()),
            round(mean(df["Length"]),2),
            round(median(df["Length"]),2),
            round(df["GC(%)"].mean(),2),
            int(df["N_Count"].sum()),
            len(invalid),
            len(dup)
        ]
    })

    with pd.ExcelWriter(out/"QC_report.xlsx",engine="openpyxl") as w:
        summary.to_excel(w,sheet_name="Summary",index=False)
        df.to_excel(w,sheet_name="Gene Statistics",index=False)
        dup.to_excel(w,sheet_name="Duplicate Sequences",index=False)
        invalid.to_excel(w,sheet_name="Invalid Sequences",index=False)

    df.to_csv(out/"gene_statistics.csv",index=False)
    dup.to_csv(out/"duplicate_sequences.csv",index=False)
    invalid.to_csv(out/"invalid_sequences.csv",index=False)

    plt.figure(figsize=(8,5))
    plt.hist(df["Length"],bins=50)
    plt.title("Gene Length Distribution")
    plt.xlabel("Length (bp)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out/"length_distribution.png")
    plt.close()

    plt.figure(figsize=(8,5))
    plt.hist(df["GC(%)"],bins=50)
    plt.title("GC Distribution")
    plt.xlabel("GC (%)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out/"gc_distribution.png")
    plt.close()

    with open(out/"QC_summary.txt","w") as f:
        f.write(summary.to_string(index=False))

    print(summary)
    print("\nQC complete.")
    print(f"Results saved to: {out.resolve()}")

if __name__=="__main__":
    main()
