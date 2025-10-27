import ROOT
import re

def compare_and_modify_histograms(file1, file2, file3, output_root_file, modified_output_root_file):
    # X-axis labels based on histogram keywords
    x_axis_labels = {
        "pt": "p_{t} [GeV]",
        "eta": "#eta",
        "phi": "#phi [rad]",
        "E": "E [GeV]",
        "pz": " p_{z} [GeV]",
        "m": " m [GeV]",
        "Et": "Rest Energy [GeV]"
    }

    particle_labels = {
        "higgs": "Higgs Boson (h1)",
        "xd": "Dark Matter (xd)",
        "xd_": "Anti-Dark Matter (xd~)",
        "xd_xd": "DM-(Anti DM)",
        "ps_a": "Pseudoscalar (a)",
        "ps_A": "Heavy Pseudoscalar (A)",
        "bj": "Bottom Quark Jet",
        "bj_": "Bottom Anti_Quark Jet",
    }

    extra_params = [
        "M_{a} = 500 GeV",
        "M_{#chi}= 10 GeV"
    ]

    # Step 1: Compare histograms and save canvases into ROOT
    def compare_histograms(file1, file2, file3, output_root_file):
        f1 = ROOT.TFile(file1)
        f2 = ROOT.TFile(file2)
        f3 = ROOT.TFile(file3)
        output_file = ROOT.TFile(output_root_file, "RECREATE")

        # Automatically get histogram names from the first file
        hist_names = [key.GetName() for key in f1.GetListOfKeys()
                      if isinstance(f1.Get(key.GetName()), ROOT.TH1)]

        for hist_name in hist_names:
            hist1 = f1.Get(hist_name)
            hist2 = f2.Get(hist_name)
            hist3 = f3.Get(hist_name)

            if not (hist1 and hist2 and hist3):
                print(f"Skipping '{hist_name}': not found in all files.")
                continue

            # --- normalization ---#
            if hist1.Integral() != 0:
                hist1.Scale(0.03401011318)
            if hist2.Integral() != 0:
                hist2.Scale(0.7855705028000001)
            if hist3.Integral() != 0:
                hist3.Scale(0.7213413323000001)
            #--------------

            # Create square canvas and styling
            canvas = ROOT.TCanvas(hist_name, hist_name, 800, 800)
            # Set margins suited for square plotting frame
            canvas.SetLeftMargin(0.15)
            canvas.SetRightMargin(0.15)
            canvas.SetTopMargin(0.10)
            canvas.SetBottomMargin(0.15)
            canvas.SetTicks(1, 1)
            canvas.SetLogy()  # keep log scale

            # Prepare histogram styles (use hist2 as frame/reference)
            hist2.SetLineColor(ROOT.kBlue)
            hist1.SetLineColor(ROOT.kRed)
            hist3.SetLineColor(ROOT.kGreen)

            hist1.SetStats(0)
            hist2.SetStats(0)
            hist3.SetStats(0)

            # Make sure hist2 provides the axes/frame so we can control axis ranges and titles
            hist2.GetYaxis().SetTitle("Number of Events")
            hist2.GetYaxis().CenterTitle(True)
            hist2.GetXaxis().CenterTitle(True)

            # Draw frame first and then others
            hist2.Draw("HIST")
            hist1.Draw("HIST SAME")
            hist3.Draw("HIST SAME")

            # Set X axis title if keyword matches
            for key, label in x_axis_labels.items():
                if key in hist_name:
                    hist2.GetXaxis().SetTitle(label)
                    break

            # Manual Y-range padding (safe for log scale)
            # Get reasonable non-zero minimum
            y_min = hist2.GetMinimum()
            y_max = hist2.GetMaximum()
            # If histogram contains zeros or negative values (common with log scale),
            # choose a small positive floor for the minimum to display in log scale.
            if y_min <= 0 or ROOT.gPad.GetLogy():
                # find the smallest non-zero bin content among hist2, hist1, hist3
                min_nonzero = None
                for h in (hist1, hist2, hist3):
                    for ibin in range(1, h.GetNbinsX() + 1):
                        val = h.GetBinContent(ibin)
                        if val > 0:
                            if min_nonzero is None or val < min_nonzero:
                                min_nonzero = val
                if min_nonzero is None:
                    # fallback tiny positive value
                    min_nonzero = 1e-6
                # set lower bound a bit below min_nonzero but positive
                y_floor = max(min_nonzero * 0.5, 1e-9)
                y_top = max(y_max * 1.2, min_nonzero * 10.0)
                hist2.SetMinimum(y_floor)
                hist2.SetMaximum(y_top)
            else:
                # linear scale safe padding
                hist2.SetMinimum(y_min * 0.9 if y_min > 0 else y_min - abs(y_min) * 0.1)
                hist2.SetMaximum(y_max * 1.2)

            # Manual X-range padding similar to Step 2
            bin_min = hist2.GetXaxis().GetFirst()
            bin_max = hist2.GetXaxis().GetLast()
            x_min = hist2.GetBinLowEdge(bin_min)
            x_max = hist2.GetXaxis().GetBinUpEdge(bin_max)
            # Protect against degenerate axis
            if x_max == x_min:
                # fallback to histogram's actual axis min/max
                x_min = hist2.GetXaxis().GetXmin()
                x_max = hist2.GetXaxis().GetXmax()
            padding = 0.01 * (x_max - x_min) if (x_max - x_min) != 0 else 0.01 * abs(x_min if x_min != 0 else 1.0)
            hist2.GetXaxis().SetRangeUser(x_min - padding, x_max + padding)
            
            # Particle + extra params label using TLatex (top-left **outside** frame, in ONE row)
            latex = ROOT.TLatex()
            latex.SetNDC(True)
            latex.SetTextFont(43)
            latex.SetTextSize(20)
            latex.SetTextAlign(13)

# Start near the top-left, outside plotting frame
            x_pos = 0.16
            y_pos = 0.94

# First draw particle label
            for key, label in particle_labels.items():
                if key in hist_name:
                    latex.DrawLatex(x_pos, y_pos, label)
                    x_pos += 0.29   # shift right for next text block
                    break

# Then draw all extra params on the same row
            for param in extra_params:
                latex.DrawLatex(x_pos, y_pos, param)
                x_pos += 0.20  # move further right each time

            
            
           

            # Legend
            legend = ROOT.TLegend(0.70, 0.70, 0.78, 0.85)
            legend.SetBorderSize(0)
            legend.SetTextFont(43)
            legend.SetTextSize(20)
            legend.AddEntry(hist1, "tan#theta = 1.0", "l")
            legend.AddEntry(hist2, "tan#theta = 10.0", "l")
            legend.AddEntry(hist3, "tan#theta = 15.0", "l")
            legend.Draw()

            # Enforce square plotting frame (equal aspect ratio)
            canvas.Update()  # must update before touching gPad
            try:
                # Use a fixed aspect ratio of 1:1 for the pad plotting area
                ROOT.gPad.SetFixedAspectRatio(1.0)
            except Exception:
                # If SetFixedAspectRatio not available in your ROOT, fall back to adjusting margins
                pass

            # Re-draw to make sure axes apply
            canvas.Modified()
            canvas.Update()

            canvas.Write()
            print(f"Saved comparison for histogram '{hist_name}'.")

        f1.Close()
        f2.Close()
        f3.Close()
        output_file.Close()

    # Step 2: Modify histograms + save to new ROOT + multipage PDF
    def modify_histograms(input_root_file, output_root_file):
        input_file = ROOT.TFile(input_root_file, "READ")
        output_file = ROOT.TFile(output_root_file, "RECREATE")
        keys = input_file.GetListOfKeys()

        c = ROOT.TCanvas("c", "c", 800, 800)
        c.SetLeftMargin(0.15)
        c.SetRightMargin(0.15)
        c.SetTopMargin(0.1)
        c.SetBottomMargin(0.15)
        c.SetTicks(1, 1)

        c.Print("output_tan_norm.pdf[")

        for key in keys:
            obj_name = key.GetName()
            canvas = input_file.Get(obj_name)

            if not isinstance(canvas, ROOT.TCanvas):
                continue

            primitives = canvas.GetListOfPrimitives()
            for primitive in primitives:
                # Match TH1 and TH2 improvements; you used TH2 earlier but be flexible
                if isinstance(primitive, ROOT.TH2) or isinstance(primitive, ROOT.TH1):
                    # set y-axis title for consistency
                    try:
                        primitive.GetYaxis().SetTitle("Number of Events")
                        primitive.GetYaxis().CenterTitle(True)
                        primitive.GetXaxis().CenterTitle(True)
                    except Exception:
                        pass

                    for key_, label in x_axis_labels.items():
                        if key_ in primitive.GetName():
                            primitive.GetXaxis().SetTitle(label)
                            break

                    canvas.SetCanvasSize(800, 800)
                    canvas.SetWindowSize(900, 900)
                    primitive.GetXaxis().SetTitleOffset(1.3)
                    primitive.GetYaxis().SetTitleOffset(1.5)
                    primitive.GetXaxis().SetLabelSize(0.04)
                    primitive.GetYaxis().SetLabelSize(0.04)
                    primitive.GetXaxis().SetTitleSize(0.045)
                    primitive.GetYaxis().SetTitleSize(0.045)
                    primitive.SetLineWidth(3)

                    # Get histogram min and max from non-empty bins
                    try:
                        y_min = primitive.GetMinimum()
                        y_max = primitive.GetMaximum()
                        if y_min <= 0 or ROOT.gPad.GetLogy():
                            # compute smallest positive bin content
                            min_nonzero = None
                            nb = primitive.GetNbinsX()
                            for ibin in range(1, nb + 1):
                                val = primitive.GetBinContent(ibin)
                                if val > 0:
                                    if min_nonzero is None or val < min_nonzero:
                                        min_nonzero = val
                            if min_nonzero is None:
                                min_nonzero = 1e-6
                            primitive.SetMinimum(max(min_nonzero * 0.5, 1e-9))
                            primitive.SetMaximum(max(y_max * 1.2, min_nonzero * 10.0))
                        else:
                            primitive.SetMinimum(y_min * 0.9 if y_min > 0 else y_min - abs(y_min) * 0.1)
                            primitive.SetMaximum(y_max * 1.2)
                    except Exception:
                        pass

                    # X-axis padding
                    try:
                        bin_min = primitive.GetXaxis().GetFirst()
                        bin_max = primitive.GetXaxis().GetLast()
                        x_min = primitive.GetBinLowEdge(bin_min)
                        x_max = primitive.GetXaxis().GetBinUpEdge(bin_max)
                        padding = 0.01 * (x_max - x_min)
                        primitive.GetXaxis().SetRangeUser(x_min - padding, x_max + padding)
                    except Exception:
                        pass

                    primitive.SetTitle("g g > h1 xd xd~")

            canvas.Update()
            canvas.Write()
            canvas.Print("output_tan_norm.pdf")  # add one page to PDF

        # Close multipage PDF
        c.Print("output_tan_norm.pdf]")

        input_file.Close()
        output_file.Close()
        print(f"Modified histograms saved to '{output_root_file}' and to 'output_tan_norm.pdf'.")

    # Run both steps
    compare_histograms(file1, file2, file3, output_root_file)
    modify_histograms(output_root_file, modified_output_root_file)


# Example usage (unchanged)
if __name__ == "__main__":
    compare_and_modify_histograms(
        "../bin/template_gg_h1xdxd_sin/Events/run_06_decayed_1/unweighted_events.root",
        "../bin/template_gg_h1xdxd_sin/Events/run_14_decayed_1/unweighted_events.root",
        "../bin/template_gg_h1xdxd_sin/Events/run_15_decayed_1/unweighted_events.root",
        "histogram_tan_comparison.root",
        "modified_histogram_tan_comparison.root"
    )

