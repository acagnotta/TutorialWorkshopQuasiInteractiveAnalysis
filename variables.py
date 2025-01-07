import math
import numpy as np
class variable(object):
    def __init__(self, name, title, taglio=None, nbins=None, xmin=None, xmax=None, xarray=None, MConly = False, noUnOvFlowbin = False):
        self._name = name
        self._title = title
        self._taglio = taglio
        self._nbins = nbins
        self._xmin = xmin
        self._xmax = xmax
        self._xarray = xarray
        self._MConly = MConly
    def __str__(self):
        return  '\"'+str(self._name)+'\",\"'+str(self._title)+'\",\"'+str(self._taglio)+'\",'+str(self._nbins)+','+str(self._xmin)+','+str(self._xmax)

vars = []

vars.append(variable(name = "PuppiMET_T1_pt_nominal", title= "p_{T}^{miss} [GeV]", nbins = 12, xmin = 250, xmax=850))
vars.append(variable(name = "PuppiMET_T1_phi_nominal", title= "Puppi MET #phi", nbins = 6, xmin = -math.pi, xmax=math.pi))

# vars.append(variable(name = "LeadingFatJetPt_pt", title= "Leading FatJet p_{T} [GeV]", nbins = 8, xmin = 50, xmax=850))

vars.append(variable(name = "nJet", title= "# Jet", nbins = 10, xmin = -0.5, xmax=9.5))
vars.append(variable(name = "nFatJet", title= "# FatJet", nbins = 5, xmin = -0.5, xmax=4.5))
vars.append(variable(name = "MinDelta_phi", title= "min #Delta #phi", nbins = 18, xmin = 0, xmax = math.pi))
vars.append(variable(name = "PV_npvsGood", title= "Number of PV", nbins = 25, xmin = -0.5, xmax = 49.5))
# vars.append(variable(name = "Top_score", title= "Top Score", nbins = 40, xmin = 0, xmax=1))

vars.append(variable(name = "Top_mass", title= "Top mass [GeV]", nbins = 30, xmin = 100, xmax=250))
vars.append(variable(name = "Top_pt", title= "Top p_{T} [GeV]", nbins = 30, xmin = 100, xmax=1000))
vars.append(variable(name = "Top_score", title= "Top Score", nbins = 40, xmin = 0, xmax=1))
vars.append(variable(name = "MT_T", title= "M_{T} [GeV]", nbins = 30, xmin = 500, xmax=2000))


regions = {

    "Presel"           : "",
    "SR"               : "Top_score>0.6",
    "SR0fjets"         : "Top_score>0.6 && nForwardJet==0",
    "SRatleast1fjets"  : "Top_score>0.6 && nForwardJet>0",
    
}