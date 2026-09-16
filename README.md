# Differential Expression Analyzer GUI 🧬🖥️

A user-friendly graphical interface built with Python for performing **differential gene expression (DGE) analysis** from RNA-seq data. It supports both **TPM/normalized** and **raw count** input formats and provides visualizations like volcano plots and heatmaps.

---

## ✨ Features

- 📁 Upload gene expression files in `.csv`, `.tsv`, or `.txt` formats  
- 🔍 Specify two groups of samples for DGE comparison  
- ⚙️ Supports both TPM/normalized and raw count data  
- 📊 Visualizes results with **volcano plot** and **heatmap of top DEGs**  
- 🧠 Automatically filters low-expression genes  
- 📤 Export significant DEGs to CSV  
- 🎨 Interactive matplotlib plots inside the GUI

---

## 📦 Requirements

Install required packages using pip:

```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn scikit-learn
