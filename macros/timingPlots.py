# Macro moved to SWAN for quick studies
# SWAN_projects/CLUE/timing.ipynb

from __future__ import print_function
from ROOT import TCanvas, TGraph, TH1F, gROOT
import ROOT

import tdrstyle
tdrstyle.setTDRStyle()

from array import array
import pandas 
import pdb
import re

CLUEALGO_NAME = "ClueGaudiAlgorithmWrapper"
TOPOALGO_NAME = "TopoClusters"
SW_NAME = "CreateCaloClustersSlidingWindow"

def readFile(fileName):
    f = open(fileName, "r")
    print(fileName)

    clue_nums_hits = array( 'd' )
    clue_times = array( 'd' )
    topo_nums_hits = array( 'd' )
    topo_times = array( 'd' )
    sw_nums_hits = array( 'd' )
    sw_times = array( 'd' )

    lines = 0
    tileSizeSearch = True
    ts_tot = 0
    ts_cols = 0
    ts_rows = 0
    ts_time = 0
    for line in f:
        if "tiles size (cols,rows):" in line and tileSizeSearch:
            tileSizeSearch = False
            tileSize = re.findall('[0-9]+', line)
            ts_tot = int(tileSize[0])
            ts_cols = int(tileSize[1])
            ts_rows = int(tileSize[2])
        if "Set up time (Barrel):" in line:
            tileTime = re.findall('[0-9.]+', line)
            ts_time = float(tileTime[0])

        if "Number of calo hits:" in line:
            num_hits = re.findall('[0-9]+', line)
            if CLUEALGO_NAME in line:
                clue_nums_hits.append(int(num_hits[0]))
            elif TOPOALGO_NAME in line:
                topo_nums_hits.append(int(num_hits[0]))
                sw_nums_hits.append(int(num_hits[0]))

        elif "Elapsed time:" in line:
            time = re.findall('[0-9.]+', line)
            if CLUEALGO_NAME in line:
                clue_times.append(float(time[0]))
            elif TOPOALGO_NAME in line:
                topo_times.append(float(time[0]))
            elif SW_NAME in line:
                sw_times.append(float(time[0]))


        lines = lines+1

    print("Lines in file:",lines)
    f.close()

    data = ("CLUE (GaudiWrapper)", clue_nums_hits, clue_times,
          "TopoClustering", topo_nums_hits, topo_times,
          "SlidingWindowClustering", sw_nums_hits, sw_times,
          "Tile Size (tot, col, row)", ts_tot, ts_cols, ts_rows, ts_time)

    if len(data[1]) != len(data[2]):
        print("ERROR: data in input have different sizes!")

    return data

def plotTimingHisto(data):

    c = TCanvas('c', 'Hits vs time' )
    c.SetFillColor(0)
    minX = 0
    maxX = 30
    nBins = 120
    h1 = TH1F("","", nBins, minX, maxX)
    h1.SetLineWidth(2)
    h1.SetFillColor(ROOT.kRed+1)
    h1.SetLineColor(ROOT.kRed+1)
    h1.SetFillStyle(3004)
    for time in data[2]:
        h1.Fill(time)
    h1.GetYaxis().SetTitle( 'a.u.' )
    h1.GetXaxis().SetTitle( 'Clustering Time [ms]' )
    h1.Scale(1./h1.Integral());
    h1.Draw( 'HIST' )
    
    h2 = TH1F("","", nBins, minX, maxX)
    h2.SetLineWidth(2)
    h2.SetLineColor(ROOT.kCyan+2)
    h2.SetFillColor(ROOT.kCyan+2)
    h2.SetFillStyle(3001)
    for time in data[8]:
        h2.Fill(time)
    h2.Scale(1./h2.Integral());
    h2.Draw( 'HISTSAME' )
    
    h3 = TH1F("","", nBins, minX, maxX)
    h3.SetLineWidth(2)
    h3.SetLineColor(ROOT.kCyan+4)
    h3.SetFillColor(ROOT.kCyan+4)
    h3.SetFillStyle(3003)
    for time in data[5]:
        h3.Fill(time)
    h3.Scale(1./h3.Integral());
    h3.Draw( 'HISTSAME' )
    
    legend_histo = ROOT.TLegend(0.5,0.6,0.85,0.85)
    header = "LAr, 10 GeV gamma";
    legend_histo.SetHeader(header)
    legend_histo.AddEntry(h1, CLUEALGO_LABEL, 'f')
    legend_histo.AddEntry(h2, SWALGO_LABEL, 'f')
    legend_histo.AddEntry(h3, TOPOALGO_LABEL, 'f')
    legend_histo.Draw('same')
    
    c.Draw()
    c.SaveAs(FOLDER+"/plots_timingStudies16/timeComparison.png", "png")
    c.SaveAs(FOLDER+"/plots_timingStudies16/timeComparison.pdf", "pdf")
    c.SaveAs(FOLDER+"/plots_timingStudies16/timeComparison.eps", "eps")

    return 0

def plotTimingScatterPlot(data):
    c = TCanvas('c', 'Hits vs time', 200, 10, 700, 500 )
    
    c.SetFillColor(0)
    c.SetGrid()
    #c.SetLogx()
    #c.SetLogy()
    
    gr = TGraph( len(data[1]), data[2], data[1] )
    gr.GetXaxis().SetLimits(0,60)
    gr.SetMarkerColor(ROOT.kRed+1)
    gr.SetMarkerStyle( 21 )
    gr.GetYaxis().SetTitle( '#Hits' )
    gr.GetXaxis().SetTitle( 'Clustering Time [ms]' )
    gr.Draw( 'AP' )
    
    gr2 = TGraph( len(data[7]), data[8], data[7] )
    gr2.SetMarkerColor(ROOT.kCyan+2)
    gr2.SetMarkerStyle( 25 )
    gr2.Draw( 'psame' )
    
    gr3 = TGraph( len(data[4]), data[5], data[4] )
    gr3.SetMarkerColor(ROOT.kCyan+4)
    gr3.SetMarkerStyle( 22 )
    gr3.Draw( 'psame' )
    
    legend_graph = ROOT.TLegend(0.5,0.2,0.89,0.45)
    header = "LAr, 10 GeV gamma";
    legend_graph.SetHeader(header)
    legend_graph.AddEntry(gr, CLUEALGO_LABEL, 'P')
    legend_graph.AddEntry(gr2, SWALGO_LABEL, 'P')
    legend_graph.AddEntry(gr3, TOPOALGO_LABEL, 'P')
    legend_graph.Draw('same')
    
    c.Draw()
    c.SaveAs(FOLDER+"/plots_timingStudies16/hitsVsTimeComparison.png", "png")
    c.SaveAs(FOLDER+"/plots_timingStudies16/hitsVsTimeComparison.pdf", "pdf")
    c.SaveAs(FOLDER+"/plots_timingStudies16/hitsVsTimeComparison.eps", "eps")
    return 0

if __name__ == "__main__":
    CLUEALGO_LABEL = "CLUE Clusters"
    TOPOALGO_LABEL = "Topo Clusers"
    SWALGO_LABEL = "Sliding Window Clusters"
    FOLDER = "/eos/home-e/ebrondol/SWAN_projects/CLUE/data/lar/timingStudies/"
    fileName = FOLDER + "log_output_gamma_10GeV_barrel_timingStudies16.log"
    data = readFile(fileName)
    plotTimingHisto(data)
    plotTimingScatterPlot(data)

