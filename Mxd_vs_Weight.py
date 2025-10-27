import ROOT
import sys
import os
import re

# Parameters (optional)
mh2 = mh3 = mhc = 700

# Read filename from command line or default
filename = sys.argv[1] if len(sys.argv) > 1 else "x_section"

# Automatically add .txt if missing
if not os.path.isfile(filename):
    if os.path.isfile(filename + ".txt"):
        filename += ".txt"
    else:
        print(f"❌ File not found: {filename} or {filename}.txt")
        sys.exit(1)

# Read and parse data
data = []
with open(filename, "r") as file:
    for line in file:
        line = line.strip()

        # Skip separators, headers, blank lines
        if (
            not line
            or line.startswith("+")
            or "Mxd" in line
            or "Weight" in line
            or "=" in line
        ):
            continue

        # Remove table borders like '|'
        line = line.replace("|", " ")
        # Replace multiple spaces or tabs with single space
        line = re.sub(r"\s+", " ", line)
        parts = line.split()

        # Expect two numbers: Mxd, Weight
        if len(parts) >= 2:
            try:
                mxd = float(parts[0])
                weight = float(parts[1])
                data.append((mxd, weight))
            except ValueError:
                continue

# Check that we got data
if not data:
    print(f"❌ No valid numeric data found in {filename}")
    sys.exit(1)

# Create TGraph
graph = ROOT.TGraph(len(data))
for i, (mxd, weight) in enumerate(data):
    graph.SetPoint(i, mxd, weight)

# Style
graph.SetTitle(f"Cross section vs Mxd (mh2=mh3=mhc={mh2} GeV);Mxd;Weight [pb]")
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1)
graph.SetMarkerColor(ROOT.kBlue)
graph.SetLineColor(ROOT.kBlue)

# Draw
c = ROOT.TCanvas("c", "Mxd vs Weight", 800, 600)
graph.Draw("APL")

# Save plot
c.SaveAs("Mxd_vs_Weight.png")
c.SaveAs("Mxd_vs_Weight.pdf")

print(f"✅ Plot saved as Mxd_vs_Weight.png and Mxd_vs_Weight.pdf using data from {filename}")

