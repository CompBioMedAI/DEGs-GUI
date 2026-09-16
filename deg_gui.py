#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os

class DEGApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Differential Expression Analyzer")
        self.root.geometry("500x500")
        self.root.configure(bg="#e6f2ff")  # Changed background color

        self.file_path = None
        self.df = None

        tk.Label(root, text="Differential Expression Analyzer", font=("Arial", 16, "bold"), bg="#e6f2ff").pack(pady=15)
        tk.Button(root, text="Upload Gene Expression File", font=("Arial", 10), command=self.upload_file).pack(pady=10)

        tk.Label(root, text="Data Type:", bg="#e6f2ff").pack()
        self.data_type_var = tk.StringVar(value="TPM")
        tk.Radiobutton(root, text="TPM / Normalized", variable=self.data_type_var, value="TPM", bg="#e6f2ff").pack()
        tk.Radiobutton(root, text="Raw Count (use DESeq2/edgeR)", variable=self.data_type_var, value="Counts", bg="#e6f2ff").pack()

        tk.Label(root, text="Group 1 Samples (comma-separated):", bg="#e6f2ff").pack()
        self.group1_entry = tk.Entry(root, width=50)
        self.group1_entry.pack(pady=5)

        tk.Label(root, text="Group 2 Samples (comma-separated):", bg="#e6f2ff").pack()
        self.group2_entry = tk.Entry(root, width=50)
        self.group2_entry.pack(pady=5)

        tk.Label(root, text="P-value Threshold:", bg="#e6f2ff").pack()
        self.pvalue_entry = tk.Entry(root, width=10)
        self.pvalue_entry.insert(0, "0.05")
        self.pvalue_entry.pack(pady=5)

        tk.Button(root, text="Find DEGs", font=("Arial", 10), bg="#4CAF50", fg="white", command=self.find_degs).pack(pady=20)

    def upload_file(self):
        self.file_path = filedialog.askopenfilename(
            initialdir=os.path.expanduser("~"),
            title="Select Gene Expression File",
            filetypes=[
                ("All Supported", "*.csv *.CSV *.tsv *.TSV *.txt *.TXT"),
                ("CSV files", "*.csv *.CSV"),
                ("TSV files", "*.tsv *.TSV"),
                ("Text files", "*.txt *.TXT"),
                ("All files", "*.*")
            ]
        )

        if self.file_path:
            try:
                if self.file_path.endswith(('.tsv', '.txt', '.TSV', '.TXT')):
                    self.df = pd.read_csv(self.file_path, sep='\t')
                else:
                    self.df = pd.read_csv(self.file_path)

                first_col = self.df.columns[0]
                if first_col != "Gene":
                    self.df.rename(columns={first_col: "Gene"}, inplace=True)

                if "Gene" not in self.df.columns:
                    messagebox.showerror("Error", "First column must contain gene names.")
                    return

                tk.Label(self.root, text=f"Loaded: {os.path.basename(self.file_path)}", bg="#e6f2ff", font=("Arial", 9, "italic")).pack()
                messagebox.showinfo("Success", f"File '{os.path.basename(self.file_path)}' loaded successfully.")

            except Exception as e:
                messagebox.showerror("File Error", str(e))

    def find_degs(self):
        if self.df is None:
            messagebox.showerror("Error", "Please upload a dataset first.")
            return

        try:
            group1 = [col.strip() for col in self.group1_entry.get().split(",")]
            group2 = [col.strip() for col in self.group2_entry.get().split(",")]
            p_thresh = float(self.pvalue_entry.get())
            data_type = self.data_type_var.get()

            gene_names = []
            p_values = []
            log2fcs = []

            if data_type == "TPM":
                self.df = self.df[(self.df[group1 + group2].mean(axis=1)) > 1]
                for idx, row in self.df.iterrows():
                    try:
                        g1_vals = row[group1].astype(float)
                        g2_vals = row[group2].astype(float)
                        stat, pval = ttest_ind(g1_vals, g2_vals, equal_var=False)

                        g1_mean = np.mean(g1_vals)
                        g2_mean = np.mean(g2_vals)
                        log2fc = np.log2(g2_mean / g1_mean) if g1_mean > 0 and g2_mean > 0 else np.nan

                        gene_names.append(row['Gene'])
                        p_values.append(pval)
                        log2fcs.append(log2fc)
                    except Exception:
                        continue

                adj_pvals = multipletests(p_values, alpha=p_thresh, method='fdr_bh')[1]

                degs = pd.DataFrame({
                    "Gene": gene_names,
                    "log2FoldChange": log2fcs,
                    "P-Value": p_values,
                    "Adjusted P-Value": adj_pvals
                })

                degs_filtered = degs[degs["Adjusted P-Value"] <= p_thresh]

                out_path = os.path.splitext(self.file_path)[0] + "_DEGs.csv"
                degs_filtered.to_csv(out_path, index=False)
                messagebox.showinfo("DEGs Saved", f"{len(degs_filtered)} significant DEGs saved to {out_path}")

            else:
                messagebox.showinfo("Raw Count Notice", "Raw count-based DE analysis (DESeq2/edgeR) is not yet implemented in this GUI. Please use raw count files in external tools.")

        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = DEGApp(root)
    root.mainloop()


