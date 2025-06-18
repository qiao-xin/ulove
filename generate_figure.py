#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Author: Xin Qiao
# Date: 27 Mar 2025

import matplotlib.pyplot as plt
import pandas as pd
import re
from sys import argv

try:
    argv[1]
    argv[2]
    argv[3]
except:
    print ("Usage: python %s <species code> <ulove summary file> <output directory>" % argv[0])
    sys.exit()

# Function to parse ULOVE summary files
def parse_busco_file(file_path):
    with open(file_path, "r") as f:
        for line in f:
            match = re.search(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)", line)
            if match:
                return {
                    "Complete (C)": float(match.group(1)),
                    "Complete and single-copy ULOVEs (S)": float(match.group(2)),
                    "Complete and duplicated ULOVEs (D)": float(match.group(3)),
                    "Fragmented ULOVEs (F)": float(match.group(4)),
                    "Missing ULOVEs (M)": float(match.group(5)),
                    "Total ULOVE groups searched": int(match.group(6)),
                }

# List of species and corresponding ULOVE files
busco_files = {
    argv[1]: argv[2],
}

# Read data from files
data = {species: parse_busco_file(file) for species, file in busco_files.items()}

# Convert to DataFrame
df = pd.DataFrame.from_dict(data, orient="index")

# Drop "Complete (C)" column as requested
df_filtered = df.drop(columns=["Complete (C)", "Total ULOVE groups searched"])

# Define colors for the categories
colors = {
    "Complete and single-copy ULOVEs (S)": "#56B4E9",  # Light Blue
    "Complete and duplicated ULOVEs (D)": "#3492C7",  # Darker Blue
    "Fragmented ULOVEs (F)": "#F0E442",  # Yellow
    "Missing ULOVEs (M)": "#F04442",  # Red
}

# Function to format ULOVE result text
def format_busco_text(row):
    return f"C:{row['Complete (C)']}%[S:{row['Complete and single-copy ULOVEs (S)']}%,D:{row['Complete and duplicated ULOVEs (D)']}%],F:{row['Fragmented ULOVEs (F)']}%,M:{row['Missing ULOVEs (M)']}%,n:{int(row['Total ULOVE groups searched'])}"

# Create formatted text for each species
df["ULOVE_Text"] = df.apply(format_busco_text, axis=1)

# Plot the ULOVE assessment results with text annotations
fig, ax = plt.subplots(figsize=(8, 1))
df_filtered.plot(kind="barh", stacked=True, color=[colors[col] for col in df_filtered.columns], ax=ax, width=0.5)

# Remove top and right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Labels and title with adjusted spacing
ax.set_xlabel("% ULOVEs", fontname="Arial", fontsize=10)
ax.set_ylabel("")  # Remove y-label
ax.set_title("ULOVE Assessment Results", pad=30, fontname="Arial", fontsize=10)  # Adjust title padding

# Adjust legend (no title, top position, two columns)
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.45),
    ncol=2,
    frameon=False,
    prop={"family": "Arial", "size": 7},
)

# Remove grid lines
ax.grid(False)

# Ensure all bars are the same length by setting x-axis limit to 100%
ax.set_xlim(0, 100)

# Set external tick directions
ax.tick_params(axis="x", direction="out")
#ax.tick_params(axis="y", direction="out", length=0)
#plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), fontsize=10, fontname="Arial")

# Add flush-left text annotations inside the bars
for index, (species, row) in enumerate(df.iterrows()):
    ax.text(
        2,  # Small margin from the left inside the bar
        index, 
        row["ULOVE_Text"], 
        ha="left", va="center",  # Align text flush left
        fontname="Arial", fontsize=8, color="black"
    )

# Adjust layout for better fit
#plt.tight_layout()

# Save as PDF
plt.savefig("%s/%s_ulove_assessment_result.pdf" % (argv[3], argv[1]), format="pdf", bbox_inches="tight")