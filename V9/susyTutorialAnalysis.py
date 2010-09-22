#! /usr/bin/env python

import ROOT
import sys
from DataFormats.FWLite import Events, Handle

files = ["SUSYPAT.root"]
events = Events (files)
handlePF  = Handle ("std::vector<pat::MET>")
handleTC  = Handle ("std::vector<pat::MET>")
handleCalo  = Handle ("std::vector<pat::MET>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
labelPF = ("patMETsPF")
labelTC = ("patMETsTC")
labelCalo = ("patMETsAK5Calo")

f = ROOT.TFile("susy_plot.root", "RECREATE")
f.cd()

#book histogram for all MET types
metPF  = ROOT.TH1F("metPF", "pfMET",    35,  0.,350.)
metTC  = ROOT.TH1F("metTC", "tcMET",    35,  0.,350.)
metCalo  = ROOT.TH1F("metCalo", "caloMET",    35,  0.,350.)

# loop over events
i = 0
for event in events:
    i = i + 1
    print  i
    # use getByLabel, just like in cmsRun
    event.getByLabel (labelPF, handlePF)
    event.getByLabel (labelTC, handleTC)
    event.getByLabel (labelCalo, handleCalo)
    # get the products
    pfMET = handlePF.product()
    tcMET = handleTC.product()
    caloMET = handleCalo.product()
    #
    metPF.Fill( pfMET[0].pt() )
    metTC.Fill( tcMET[0].pt() )
    metCalo.Fill( caloMET[0].pt() )


f.cd()

metPF.Write()
metTC.Write()
metCalo.Write()

#Make a overlaid plot of all met collections
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetOptStat(0)
canvas = ROOT.TCanvas( 'canvas', 'Canvas with met collections', 600, 600 )
canvas.SetLogy()
metPF.SetLineColor(ROOT.kRed)
metTC.SetLineColor(ROOT.kGreen)
metCalo.SetLineColor(ROOT.kBlue)
metPF.GetYaxis().SetTitle("Entries")
metPF.GetXaxis().SetTitle("MET [GeV]")
metPF.Draw("HIST")
metTC.Draw("HISTSAME")
metCalo.Draw("HISTSAME")
l = ROOT.TLegend(0.7, 0.75, 0.92, 0.92)
l.AddEntry(metPF,"pfMET","l")
l.AddEntry(metTC,"tcMET","l")
l.AddEntry(metCalo,"caloMET","l")
l.Draw()
canvas.Print("susy_met_plot.png")


f.Close()

