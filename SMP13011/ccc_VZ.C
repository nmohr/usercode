#include "TStyle.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TBox.h"

void ccc_VZ(const char *poi = "r", double rMax=17, const char *filename="ccc_VZ.pdf") {
  //   setTDRStyle();
    gStyle->SetPadTickX(0);
    //gStyle->SetPadTickY(0);
    gStyle->SetTickLength(0., "Y");

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.13);
    c1->SetGridx(0);


   int  nChann =3;
    TH2F frame("analysis topologies",";#scale[0.8]{Best fit #mu}",1,0,2.1,nChann+2,0,nChann+2);
    frame.GetYaxis()->SetLabelOffset(-0.08);
    frame.GetXaxis()->SetTitleOffset(0.8);
   


    TGraphAsymmErrors points(nChann);
   

    for (int iChann=0; iChann < nChann; iChann++){
   
      if (iChann==2) {

      TString channel = "XXX";         
      
       channel.ReplaceAll("XXX",TString::Format("#splitline{#scale[0.95]{0-lepton}  }{#scale[0.7]{#mu = %.2f #pm %.2f}}",0.9412,0.2618));   
     
      points.SetPoint(iChann,  0.9412, iChann+0.5);
      points.SetPointError(iChann, 0.2359 , +0.2618, 0, 0);                                                                                 
	
      // c1->SetTickx(0);
  frame.GetYaxis()->SetBinLabel(iChann+1, channel);  
      }

      if (iChann==0) {

      TString channel = "XXX";         
      
      channel.ReplaceAll("XXX",TString::Format("#splitline{#scale[0.95]{2-lepton}  }{#scale[0.7]{#mu = %.2f #pm %.2f}}",1.2822,0.6132));   

      points.SetPoint(iChann,  1.2822, iChann+0.5);
      points.SetPointError(iChann, 0.5668 , 0.6132, 0, 0);                                                                                 
      frame.GetYaxis()->SetBinLabel(iChann+1, channel);  
      }

      if (iChann==1) {

      TString channel = "XXX";         
      
      channel.ReplaceAll("XXX",TString::Format("#splitline{#scale[0.95]{1-lepton}  }{#scale[0.7]{#mu =  %.2f #pm %.2f}}",1.4693,0.4333));   
      points.SetPoint(iChann,  1.4693, iChann+0.5);
      points.SetPointError(iChann, 0.3849 , 0.4333 , 0, 0);                                                                                 
      frame.GetYaxis()->SetBinLabel(iChann+1, channel);  
      }



      }



    points.SetLineColor(kRed);
    points.SetLineWidth(5);
    points.SetMarkerStyle(21);
    points.SetMarkerSize(1.7);
    frame.GetXaxis()->SetTitleSize(0.05);
    frame.GetXaxis()->SetLabelSize(0.04);
    frame.GetXaxis()->SetNdivisions(505);
    frame.GetYaxis()->SetLabelSize(0.06);
    frame.Draw(); gStyle->SetOptStat(0);
    TBox globalFitBand(1.0864-0.2109, 0, 1.0864+0.2369, nChann);
    //globalFitBand.SetFillStyle(3013);
    globalFitBand.SetFillStyle(1001);
    //globalFitBand.SetFillColor(65);
    globalFitBand.SetFillColor(kGreen);
    globalFitBand.SetLineStyle(0);
    globalFitBand.DrawClone();
    TLine globalFitLine(1.09, 0, 1.09, nChann);
    globalFitLine.SetLineWidth(4);
    //globalFitLine.SetLineColor(214);
    globalFitLine.SetLineColor(1);
    globalFitLine.DrawClone();


    TLine sepLine(-2, 2., 4, 2. );
    sepLine.SetLineWidth(1);
    //globalFitLine.SetLineColor(214);
    sepLine.SetLineColor(1);
    sepLine.SetLineStyle(kDashed);

    //sepLine.DrawClone();



    points.Draw("P SAME");
    frame.Draw("AXIS SAME");

    TLatex l;
    l.SetNDC();
    //l.SetTextAlign(12);
    //l.SetTextFont(42);
    l.SetTextSize(0.035);
    l.DrawLatex(0.17,0.83,"#sqrt{s} = 8 TeV, L = 18.9 fb^{-1}");
      l.DrawLatex(0.17,0.73,TString::Format("#scale[0.95]{Combined}#scale[1.2]{ #mu = %.2f {}_{-%.2f}^{+%.2f}}",1.0864,0.2109,0.2369));    
   
    TLatex ll;
    ll.SetNDC();
    //ll.SetTextAlign(12);
    //ll.SetTextFont(42);
    ll.SetTextSize(0.05);
    ll.DrawLatex(0.17,0.87,"CMS");
    ll.SetTextSize(0.035);
    ll.DrawLatex(0.17,0.78,"pp #rightarrow VZ; Z #rightarrow b#bar{b} ");
    ll.SetTextSize(0.05);
    ll.DrawLatex(0.8,0.8,"(a)");





    c1->SetFillStyle(4000);
	  c1->SetFrameFillStyle(1000);
	  c1->SetFrameFillColor(0);


    c1->Print(filename);
}

