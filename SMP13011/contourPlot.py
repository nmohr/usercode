import math,ROOT,array
from myutils import TdrStyles
TdrStyles.tdrStyle()

def scanLineKiller(hist, threshold=0.1):
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    for ix in range(3,nx-1):
        thisline = 0
        for iy in range(1,ny):
            points = 0
            z = hist.GetBinContent(ix,iy)
            zd = hist.GetBinContent(ix,iy+1)
            if (math.fabs(z-zd) > 1e-6):
                continue
            zL = hist.GetBinContent(ix-1,iy)
            zR = hist.GetBinContent(ix+1,iy)
            disc = math.fabs(z-zL)/(math.fabs(zL-zR)+0.1*math.fabs(zL))
            if (disc > threshold):
                if (thisline == 0):
                    thisline = 1
                    hist.SetBinContent(ix,iy, 0.5*(zL+zR))
    for iy in range(3,ny-1):
        thisline = 0
        for ix in range(1,nx):
            points = 0
            z  = hist.GetBinContent(ix,iy)
            zd = hist.GetBinContent(ix+1,iy)
            if (math.fabs(z-zd) > 1e-6):
                continue
            zL = hist.GetBinContent(ix,iy-1)
            zR = hist.GetBinContent(ix,iy+1)
            disc = math.fabs(z-zL)/(math.fabs(zL-zR)+0.1*math.fabs(zL))
            if (disc > threshold):
                thisline = 1
                hist.SetBinContent(ix,iy, 0.5*(zL+zR))
                if (ix == nx-1):
                    hist.SetBinContent(ix+1,iy, 0.5*(zL+zR))


def stepLineKiller(hist, threshold=0.1):
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    while (True):
        theline = -1
        maxdiff = threshold
        diffs = 0
        ndiffs = 0
        for ix in range(3,nx):
            points = 0
            diff = 0
            diff2 = 0
            for iy in range(1,ny-1):
                z  = hist.GetBinContent(ix,iy)
                zL = hist.GetBinContent(ix-1,iy)
                zR = hist.GetBinContent(ix+1,iy)
                if (z > 20. or zL > 20. or zR > 20.): 
                    continue
                delta = (z - 0.5*(zL+zR))
                diff += delta 
                diff2 += delta*delta 
                points += 1
            if (points == 0):
                continue
            diff = diff/points 
            diff2 = math.sqrt((diff2/points - diff*diff))
            diffs += math.fabs(diff) 
            ndiffs += 1
            if (math.fabs(diff) > maxdiff):
                theline = ix 
                maxdiff = math.fabs(diff)
        if (theline == -1): 
            break
        diffs/=ndiffs
        if (maxdiff < 3*diffs):
            break
        for iy in range(1,ny):
            zL = hist.GetBinContent(theline-1,iy)
            zR = hist.GetBinContent(theline+1,iy)
            hist.SetBinContent(theline, iy, 0.5*(zL+zR))


def spikeKiller(hist, threshold=0.2, algo=0):
    scanLineKiller(hist,0.1)
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    if (algo == 0):
        for ix in range(2,nx):
            for iy in range(2,ny):
                points = 0
                z  = hist.GetBinContent(ix,iy)
                z1 = hist.GetBinContent(ix-1,iy-1)
                z2 = hist.GetBinContent(ix-1,iy+1)
                z3 = hist.GetBinContent(ix+1,iy+1)
                z4 = hist.GetBinContent(ix+1,iy-1)
                zminD = ROOT.TMath.Min(ROOT.TMath.Min(z1,z2),ROOT.TMath.Min(z3,z4))
                zmaxD = ROOT.TMath.Max(ROOT.TMath.Max(z1,z2),ROOT.TMath.Max(z3,z4))
                z1 = hist.GetBinContent(ix,iy-1)
                z2 = hist.GetBinContent(ix,iy+1)
                z3 = hist.GetBinContent(ix+1,iy)
                z4 = hist.GetBinContent(ix-1,iy)
                zminQ = ROOT.TMath.Min(ROOT.TMath.Min(z1,z2),ROOT.TMath.Min(z3,z4))
                zmaxQ = ROOT.TMath.Max(ROOT.TMath.Max(z1,z2),ROOT.TMath.Max(z3,z4))
                disc = ROOT.TMath.Max(ROOT.TMath.Max(z-zmaxD,zminD-z), ROOT.TMath.Max(z-zmaxQ,zminQ-z))
                if (disc > threshold):
                    #print "suspicious point in %s at x = %g, y = %g: z = %g, zminQ = %g, zmaxQ = %g, zminD = %g, zmaxD = %g, discrim = %g\n"        %(hist.GetName(), hist.GetXaxis().GetBinCenter(ix), hist.GetYaxis().GetBinCenter(iy), z, zminQ, zmaxQ, zminD, zmaxD, disc)
                    hist.SetBinContent(ix,iy, 0.25*(z1+z2+z3+z4)) 
    elif (algo == 1):
        for ix in range(1,nx):
            for iy in range(1,ny):
                z  = hist.GetBinContent(ix,iy)
                z1 = hist.GetBinContent(ix,iy-1) if iy > 1 else z
                z2 = hist.GetBinContent(ix,iy+1) if iy < ny else z
                disc = math.fabs(z - 0.5*(z1+z2)) - 3*math.fabs(z1-z2) 
                if (disc > threshold):
                    hist.SetBinContent(ix,iy, 0.5*(z1+z2))
                    continue
                z1 = hist.GetBinContent(ix+1,iy) if ix > 1 else z
                z2 = hist.GetBinContent(ix-1,iy) if ix < ny else z
                disc = math.fabs(z - 0.5*(z1+z2)) - 3*math.fabs(z1-z2)
                if (disc > threshold):
                    hist.SetBinContent(ix,iy, 0.5*(z1+z2))
                    continue
    stepLineKiller(hist)


def frameTH2D(th2):

    frameValue = 1000
    if "bayes" in th2.GetName():
        frameValue = 0.0

    xw = th2.GetXaxis().GetBinWidth(0)
    yw = th2.GetYaxis().GetBinWidth(0)

    nx = th2.GetNbinsX()
    ny = th2.GetNbinsY()

    x0 = th2.GetXaxis().GetXmin()
    x1 = th2.GetXaxis().GetXmax()

    y0 = th2.GetYaxis().GetXmin()
    y1 = th2.GetYaxis().GetXmax()

    framed = ROOT.TH2D("%s framed" %th2.GetName(),"%s framed" %th2.GetTitle(),nx + 2, x0-xw, x1+xw,ny + 2, y0-yw, y1+yw)

    #copy content
    for ix in range(0,nx):
        for iy in range(0, ny):
            framed.SetBinContent(1+ix, 1+iy, th2.GetBinContent(ix,iy))
        
    
    #Frame with huge values
    nx = framed.GetNbinsX()
    ny = framed.GetNbinsY()
    for ix in range(1,nx):
        framed.SetBinContent(ix,  1, frameValue)
        framed.SetBinContent(ix, ny, frameValue)
    for iy in range(2, ny-1):
        framed.SetBinContent( 1, iy, frameValue)
        framed.SetBinContent(nx, iy, frameValue)

    return framed

def contourFromTH2(h2in, threshold, minPoints = 20):
    contours = array.array('d')
    contours.append(threshold)

    h2 = frameTH2D(h2in)

    h2.SetContour(1, contours)

    # Draw contours as filled regions, and Save points
    h2.Draw("CONT Z LIST")
    ROOT.gPad.Update() # Needed to force the plotting and retrieve the contours in TGraphs

    # Get Contours
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    contLevel = None


    ret = ROOT.TList()
    for i in range(conts.GetSize()):
        contLevel = conts[i]
        print "Contour %d has %d Graphs" %(i,contLevel.GetSize())
        for j in range(contLevel.GetSize()):
            gr1 = contLevel.At(j)
            if (gr1.GetN() > minPoints):
                ret.Add(gr1.Clone())
    return ret

def bestFit(t, x, y):
    nfind = t.Draw(y+":"+x, "deltaNLL == 0")
    if not nfind:
        gr0 = ROOT.TGraph(1)
        gr0.SetPoint(0,-999,-999)
        gr0.SetMarkerStyle(34)
        gr0.SetMarkerSize(2.0)
    else:
        gr0 = ROOT.gROOT.FindObject("Graph").Clone()
        gr0.SetMarkerStyle(34) 
        gr0.SetMarkerSize(2.0)
        if (gr0.GetN() > 1):
            gr0.Set(1)
    return gr0

def countGridPointsFromTree(t, x, xmin = -1, xmax = -1):
    if (xmin == xmax):
        xmin = t.GetMinimum(x)
        xmax = t.GetMaximum(x)
    t.Draw("%s>>h1000(1000,%10g,%10g)" %(x,xmin-1e-4,xmax+1e-4), "deltaNLL > 0")
    h1000 = ROOT.gROOT.FindObject("h1000")
    bins = 0
    for i in range(h1000.GetNbinsX()):
        if not h1000.GetBinContent(i) == 0:
            bins += 1
    del h1000
    return bins


def treeToHist2D(t, x, y, name, xmin = -1, xmax = -1, ymin = -1, ymax = -1):
    if (xmin == xmax):
        xmin = t.GetMinimum(x)
        xmax = t.GetMaximum(x)
    if (ymin == ymax):
        ymin = t.GetMinimum(y)
        ymax = t.GetMaximum(y)
    xbins = countGridPointsFromTree(t,x,xmin,xmax)
    ybins = countGridPointsFromTree(t,y,ymin,ymax)
    dx = (xmax-xmin)/(xbins-1)
    dy = (ymax-ymin)/(ybins-1)
    xmin -= 0.5*dx 
    xmax += 0.5*dx
    ymin -= 0.5*dy
    ymax += 0.5*dy
    if (math.fabs(xmin) < 1e-5):
        xmin = 0
    if (math.fabs(xmax) < 1e-5):
        xmax = 0
    #std::cout << "In making " << name << ", guessed " << xbins << " bins for " << x << " from " << xmin << " to " << xmax << std::endl
    #std::cout << "In making " << name << ", guessed " << ybins << " bins for " << y << " from " << ymin << " to " << ymax << std::endl
    t.Draw("2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)" %(y, x, name, xbins, xmin, xmax, ybins, ymin, ymax), "deltaNLL != 0", "PROF")
    prof = ROOT.gROOT.FindObject(name+"_prof")
    h2d = ROOT.TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax)
    for ix in range(1,xbins):
        for iy in range(1,ybins):
             h2d.SetBinContent(ix, iy, prof.GetBinContent(ix,iy))
    h2d.SetDirectory(0)
    return h2d


def myText(txt="CMS Preliminary",ndcX=0.6,ndcY=0.6,size=0.7):
    ROOT.gPad.Update()
    x = ndcX 
    y = ndcY 
    text = ROOT.TLatex()
    text.SetNDC()
    print '%s %s' %(txt, text.GetTextSize()*size)
    #text.SetTextColor(12)
    text.SetTextColor(ROOT.kBlack)
    text.SetTextSize(text.GetTextSize()*size)
    text.DrawLatex(x,y,txt)
    return text

def styleMultiGraph(tmg, lineColor, lineWidth, lineStyle, fillColor):
    for i in range(tmg.GetSize()):
        g = tmg.At(i)
        g.SetLineColor(lineColor)
        g.SetLineWidth(lineWidth)
        g.SetLineStyle(lineStyle)
        g.SetFillColor(fillColor)

def getShade(cont):
    shade = cont.At(0).Clone()
    shade.Set(cont.At(1).GetN()+cont.At(1).GetN())
    for i in range(cont.At(1).GetN()):
        shade.SetPoint(i,cont.At(1).GetX()[i],cont.At(1).GetY()[i])
    for i in range(cont.At(0).GetN()):
        shade.SetPoint(cont.At(1).GetN()+i,cont.At(0).GetX()[i],cont.At(0).GetY()[i])
    return shade

def doContour(fileName):
    f = ROOT.TFile.Open(fileName)
    t = f.Get('limit')
    #xsecWZ = 33.85
    xsecWZ = 22.3
    xsecZZ = 7.7
    #8.297
    theoryWZ = xsecWZ
    theoryZZ = xsecZZ
    x = '%s*r_WZHF' %xsecWZ
    y = '%s*r_ZZHF' %xsecZZ
    name = 'Test'
    th2 = treeToHist2D(t,x,y,name,0.,4.*xsecWZ,0.,3*xsecZZ)
    spikeKiller(th2)

    cont68 = contourFromTH2(th2,1.)
    cont95 = contourFromTH2(th2, 3.84)
    theBestFit = bestFit(t,x,y)
    wzError = array.array("d",[ theoryWZ*0.405*0.7454])
    zzError = array.array("d",[ theoryZZ*0.335*0.7850])
    theBestFit = ROOT.TGraphErrors(1, theBestFit.GetX(), theBestFit.GetY(), wzError, zzError)
    theBestFit.SetMarkerStyle(34)
    theBestFit.SetMarkerSize(2)
    canvas = ROOT.TCanvas(name,name, 800, 800)
    canvas.SetFillStyle(4000)
    canvas.SetFrameFillStyle(1000)
    canvas.SetFrameFillColor(0)
    h2 = ROOT.TH2F("cvcf","#kappa_{V} (scaling of vector boson couplings) #kappa_{b} (scaling of fermion couplings)",1,0,3.*xsecWZ,1,0,2.5*xsecZZ)
    #h2 = th2.Clone()
    #h2.SetNdivisions(505,'X')
    #h2.SetNdivisions(505,'Y')
    h2.Draw('AXIG')
    h2.GetXaxis().SetTitle('#sigma_{WZ} [pb]')
    h2.GetYaxis().SetTitle('#sigma_{ZZ} [pb]')
    h2.GetYaxis().SetTitleOffset(1.2)

    styleMultiGraph(cont95,ROOT.kBlack,1,ROOT.kSolid,ROOT.kBlue-9)
    styleMultiGraph(cont68,ROOT.kBlack,1,ROOT.kSolid,ROOT.kBlue-7)
    #shade95 = getShade(cont95)
    #shade95.Draw('FSAME')
    #shade68 = getShade(cont68)
    #shade68.Draw('FSAME')
    for i in range(cont95.GetSize()):
        cont95.At(i).Draw("FL SAME")
    for i in range(cont68.GetSize()):
        cont68.At(i).Draw("FL SAME")

    wz = array.array("d",[ theoryWZ])
    zz = array.array("d",[ theoryZZ])
    wzErr = array.array("d",[ theoryWZ*0.05])
    zzErr = array.array("d",[ theoryZZ*0.05])
    sm = ROOT.TEllipse(theoryWZ, theoryZZ, 0.05*theoryWZ, 0.05*theoryZZ, 0, 360)    
    sm.SetFillStyle(0)
    sm = ROOT.TGraphErrors(1,wz,zz,wzErr,zzErr)
    #sm = ROOT.TMarker(theoryWZ,theoryZZ,21)
    sm.SetLineWidth(3)
    sm.SetLineColor(ROOT.kOrange)
    sm.SetMarkerStyle(21)
    sm.SetMarkerSize(1)
    sm.SetMarkerColor(ROOT.kOrange)
    sm.Draw('SAMEP')
    
    theBestFit.Draw('SAMEP')
    legend = ROOT.TLegend(0.58,0.77,0.83,0.92) # for many entries
    legend.SetLineWidth(2)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(4000)
    legend.SetTextFont(62)
    legend.SetTextSize(0.035)
    legend.AddEntry(theBestFit, 'Best fit (68% stat.)','PEL')
    #NIKLAS
    #legend.AddEntry(cont68.At(0), '68% CL','FL')
    #legend.AddEntry(cont95.At(0), '95% CL','FL')
    legend.AddEntry(cont68.At(0), '68% CL','F')
    legend.AddEntry(cont95.At(0), '95% CL','F')
    legend.AddEntry(sm, 'MCFM NLO','PEL')
    legend.Draw()
    #t0 = myText('VH(b#bar{b})',0.72,0.60)
    t0 = myText('CMS',0.17,0.88,1.)
    t1 = myText("#sqrt{s} = 8 TeV, L = 18.9 fb^{-1}",0.17,0.83)
    t2 = myText("pp #rightarrow VZ; Z #rightarrow b#bar{b}",0.17,0.78)
    t2 = myText("Inclusive cross sections for 60 < m_{Z} < 120 GeV",0.17,0.73)
    t3 = myText("(b)",0.8,0.65,0.95)
    canvas.Print('WZ_ZZ_2d.pdf')

doContour('data/wz_zz_scan_newLumi_addBoostUncert.root')
