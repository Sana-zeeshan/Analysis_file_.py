import ROOT

def compare_and_modify_histograms(file1, file2, file3, output_root_file, modified_output_root_file):
    # X-axis labels based on histogram keywords
    x_axis_labels = {
        "pt": "Transverse Momentum (p_t) [GeV]",
        "eta": "Pseudorapidity (#eta)",
        "phi": "Azimuthal Angle (#phi) [rad]",
        "E": "Energy of Particle (E) [GeV]",
        "pz": "Longitudinal Momentum (p_z) [GeV]",
        "m": "Mass (m) [GeV]"
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
        "M_{a} = 200 GeV",
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

            canvas = ROOT.TCanvas(hist_name, hist_name, 800, 800)
            canvas.SetLeftMargin(0.15)
            hist1.SetLineColor(ROOT.kRed)
            hist2.SetLineColor(ROOT.kBlue)
            hist3.SetLineColor(ROOT.kGreen)

            hist1.SetStats(0)
            hist2.SetStats(0)
            hist3.SetStats(0)

            hist1.Draw("HIST")
            hist2.Draw("HIST SAME")
            hist3.Draw("HIST SAME")

            for key, label in x_axis_labels.items():
                if key in hist_name:
                    hist1.GetXaxis().SetTitle(label)
                    break

            legend = ROOT.TLegend(0.7, 0.7, 0.8, 0.85)
            legend.SetBorderSize(0)
            legend.SetTextFont(43)
            legend.SetTextSize(20)
            legend.AddEntry(hist1, "M_{H} = 400 GeV", "l")
            legend.AddEntry(hist2, "M_{H} = 700 GeV", "l")
            legend.AddEntry(hist3, "M_{H} = 1000 GeV", "l")
            legend.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC(True)
            latex.SetTextFont(43)
            latex.SetTextSize(20)

            y_pos = 0.82
            # start position (top-left)
            # Particle label
            for key, label in particle_labels.items():
                if key in hist_name:
                    latex.DrawLatex(0.18, y_pos, label)
                    y_pos -= 0.05  # move down for next line
                    break
            # Extra params
            for param in extra_params:
                latex.DrawLatex(0.18, y_pos, param)
                y_pos -= 0.05

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

        # Open multipage PDF
        c = ROOT.TCanvas()
        c.Print("output_h_decay1.pdf[")

        for key in keys:
            obj_name = key.GetName()
            canvas = input_file.Get(obj_name)

            if not isinstance(canvas, ROOT.TCanvas):
                continue

            primitives = canvas.GetListOfPrimitives()
            for primitive in primitives:
                if isinstance(primitive, ROOT.TH2):
                    primitive.GetYaxis().SetTitle("Number of Events")

                    for key_, label in x_axis_labels.items():
                        if key_ in primitive.GetName():
                            primitive.GetXaxis().SetTitle(label)
                            break
                    primitive.GetXaxis().CenterTitle(True)
                    primitive.GetYaxis().CenterTitle(True)
                    primitive.SetLineWidth(3)

# Get histogram min and max from non-empty bins
                    y_min = primitive.GetMinimum()
                    y_max = primitive.GetMaximum()

# Add some padding so it looks nicer
                    primitive.SetMinimum(y_min * 0.9 if y_min > 0 else y_min - abs(y_min) * 0.1)
                    primitive.SetMaximum(y_max * 1.2)
                    
                    bin_min = primitive.GetXaxis().GetFirst()
                    bin_max = primitive.GetXaxis().GetLast()

                    x_min = primitive.GetBinLowEdge(bin_min)
                    x_max = primitive.GetXaxis().GetBinUpEdge(bin_max)

# Add ~5% padding on each side
                    padding = 0.01 * (x_max - x_min)
                    primitive.GetXaxis().SetRangeUser(x_min - padding, x_max + padding)

                    primitive.SetTitle("Dark Matter Production via gluon fusion and higgs decay")

            canvas.Update()
            canvas.Write()
            canvas.Print("output_h_decay1.pdf")  # add one page to PDF

        # Close multipage PDF
        c.Print("output_h_decay1.pdf]")

        input_file.Close()
        output_file.Close()
        print(f"Modified histograms saved to '{output_root_file}' and to 'output.pdf'.")

    # Run both steps
    compare_histograms(file1, file2, file3, output_root_file)
    modify_histograms(output_root_file, modified_output_root_file)


# Example usage
compare_and_modify_histograms(
    "../bin/template_h_decay/Events/run_01_decayed_1/unweighted_events.root",
    "../bin/template_h_decay/Events/run_02_decayed_1/unweighted_events.root",
    "../bin/template_h_decay/Events/run_03_decayed_1/unweighted_events.root",
    "histogram_h_decay1_comparison.root",
    "modified_histogram_h_decay1_comparison.root"
)

