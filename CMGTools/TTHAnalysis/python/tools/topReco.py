from CMGTools.TTHAnalysis.treeReAnalyzer import *

def mt(*x0): 
    x = [ (ix.p4() if isinstance(ix,Object) else ix) for ix in x0 ]
    ht = sum([xi.Pt() for xi in x])
    pt = sum(x[1:],x[0]).Pt()
    return sqrt(max(ht*ht-pt*pt,0))
def solveWlv(lW,met,metphi):
    MW=80.4
    a = (1 - (lW.p4().Z()/lW.p4().E())**2)
    ppe    = met * lW.pt * cos(lW.phi - metphi)/lW.p4().E()
    brk    = MW**2 / (2*lW.p4().E()) + ppe
    b      = (lW.p4().Z()/lW.p4().E()) * brk
    c      = met**2 - brk**2
    delta   = b**2 - a*c
    sqdelta = sqrt(delta)    if delta    > 0 else 0
    return [ (b + s*sqdelta)/a for s in +1,-1 ]
def toplvb(bp4,lW,met):
    metx, mety = met.Px(), met.Py()
    metz = min(solveWlv(lW,met.Pt(),met.Phi()), key = lambda pz : abs(pz))
    nup4 = ROOT.TLorentzVector(metx,mety,metz,hypot(met.Pt(),metz))
    return bp4+lW.p4()+nup4
def jp4(lep):
    ret = ROOT.TLorentzVector()
    ret.SetPtEtaPhiM(lep.pt/lep.jetPtRatio,lep.eta,lep.phi,lep.mass)
    return ret
def scalWjj(wjj):
    ret = ROOT.TLorentzVector(wjj)
    sf = 80.4/wjj.M()
    ret.SetPtEtaPhiM(wjj.Pt()*sf,wjj.Eta(),wjj.Phi(),wjj.M()*sf)
    return ret
    
class GenWjj:
    def __init__(self,debug=False):
        self.debug = debug
    def listBranches(self):
        return [ ("nGenWjj",'I'),
                 ("GenWjj_ij1", "I",20,"nGenWjj"),
                 ("GenWjj_ij2", "I",20,"nGenWjj"),
                 ("GenWjj_irj1", "I",20,"nGenWjj"),
                 ("GenWjj_irj2", "I",20,"nGenWjj"),
                 ("GenWjj_iW",  "I",20,"nGenWjj"),
                 ("GenWjj_iMom","I",20,"nGenWjj"),
                 ("GenWjj_acc","I",20,"nGenWjj") ]
    def __call__(self,event):
       ## return if not on MC
        ret = {
            'nGenWjj' : 0,
            'GenWjj_ij1' : [],
            'GenWjj_ij2' : [],
            'GenWjj_irj1' : [],
            'GenWjj_irj2' : [],
            'GenWjj_iW' : [],
            'GenWjj_iMom' : [],
            'GenWjj_acc' : [],
        }
        if event.run != 1: return ret
        # make python lists as Collection does not support indexing in slices
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets = [j for j in Collection(event,"Jet","nJet")]
        gens = [j for j in Collection(event,"GenPart","nGenPart")]
        njet = len(jets)
        Ws = []
        for iW, W in enumerate(gens):
            if abs(W.pdgId) != 24: continue
            j1, j2, ij1, ij2 = None, None, -1, -1; 
            for ij1,cj1 in enumerate(gens):
                if cj1.motherIndex == iW and abs(cj1.pdgId)<6:
                    for ij2, cj2 in enumerate(gens):
                        if ij2 <= ij1: continue
                        if cj2.motherIndex == iW and abs(cj2.pdgId)<6:
                            j1 = cj1; 
                            j2 = cj2; 
                            if (j1.pt < j2.pt): 
                                (j1,j2,ij1,ij2) = (j2,j1,ij2,ij1)
                            break
                    break
            if ij1 == -1 or ij2 == -1:
                continue
            for (ij,j) in (ij1,j1), (ij2,j2):
                rj,dr = closest(j, jets)
                if dr < 0.5: j.rec = rj
                else:        j.rec = None
            W.j1 = j1; W.j2 = j2; 
            W.jj = j1.rec.p4() + j2.rec.p4() if j1.rec and j2.rec else None
            W.idx = { 'Mom':W.motherIndex, 'W':iW, 'j1':ij1, 'j2':ij2 }
            Ws.append(W) 
            ret['nGenWjj'] += 1;
            for x in "Mom W j1 j2".split():
                ret['GenWjj_i'+x].append(W.idx[x])
            ret['GenWjj_irj1'].append(jets.index(j1.rec) if j1.rec else -1)
            ret['GenWjj_irj2'].append(jets.index(j2.rec) if j2.rec else -1)
            ret['GenWjj_acc'].append( abs(W.j1.eta) < 2.4 and abs(W.j2.eta) < 2.4 and W.j1.pt > 25 and W.j2.pt > 25 )
        if not self.debug or not Ws: return ret
        for W in Ws:
                j1 = W.j1; j2 = W.j2;
                itop = W.idx['Mom']; iW = W.idx['W']; ij1 = W.idx['j1']; ij2 = W.idx['j2']
                if itop >= 0: 
                    top = gens[itop]
                    print "mom %2d       pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( itop, top.pt, top.eta, top.phi, top.pdgId)
                print "  W  %2d    pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d charge %+.0f " % ( iW, W.pt, W.eta, W.phi, W.pdgId, W.charge)
                if j1 and j2 and j1.rec and j2.rec:
                    print " |  |     jj pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                            W.jj.Pt(), W.jj.Eta(), W.jj.Phi(), W.jj.M(), deltaR(W.jj.Eta(),W.jj.Phi(),W.eta,W.phi), W.jj.Pt()/W.pt )
                for (ij,j) in (ij1,j1), (ij2,j2):
                    if not j: continue
                    if ij == ij2: print " |  |  " 
                    print " |  +- j %2d  pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d charge %+.2f " % ( ij, j.pt, j.eta, j.phi, j.pdgId, j.charge)
                    if j.rec:
                        print " |  %s    jet pT %7.2f eta %+5.2f phi %+5.2f            btag %.3f charge %+.2f qgl %+.2f dr %5.3f ptR %4.2f" % ( " " if ij == ij2 else "|",
                            j.rec.pt, j.rec.eta, j.rec.phi, max(0,j.rec.btagCSV), j.rec.charge, j.rec.qgl, deltaR(j,j.rec), j.rec.pt/j.pt )
                print ""
        if any(W for W in Ws if W.j1 and W.j1.pt > 25 and W.j2 and W.j2.pt > 25 and (W.j1.rec is None or W.j2.rec is None)):
            print "Where did the jets go?"
            for l in leps:
                print "  lepton pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d   mcMatchId %+3d mcMatchAny %+3d " % ( 
                    l.pt, l.eta, l.phi, l.pdgId, l.mcMatchId, l.mcMatchAny)
            for j in jets:
                print "  jet    pT %7.2f eta %+5.2f phi %+5.2f btag %.3f   mcMatchId %+3d mcFlavour  %+3d " % (
                    j.pt, j.eta, j.phi, max(j.btagCSV,0.), j.mcMatchId, j.mcFlavour)
            print ""
        print "-" * 100
        print ""
        return ret

class GenTopWjj:
    def __init__(self,debug=False):
        self.debug = debug
    def listBranches(self):
        return [ ("nGenTopWjj",'I'),
                 ("GenTopWjj_iGenWjj", "I",20,"nGenTopWjj"),
                 ("GenTopWjj_ib", "I",20,"nGenTopWjj"),
                 ("GenTopWjj_itop","I",20,"nGenTopWjj"),
                 ("GenTopWjj_rhad","I",20,"nGenTopWjj"),
                 ("GenTopWjj_irb", "I",20,"nGenTopWjj") ]
    def __call__(self,event):
       ## return if not on MC
        ret = {
            'nGenTopWjj' : 0,
            'GenTopWjj_iGenWjj' : [],
            'GenTopWjj_itop' : [],
            'GenTopWjj_ib' : [],
            'GenTopWjj_rhad' : [],
            'GenTopWjj_irb' : [],
        }
        if event.run != 1: return ret
        # make python lists as Collection does not support indexing in slices
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets = [j for j in Collection(event,"Jet","nJet")]
        gens = [j for j in Collection(event,"GenPart","nGenPart")]
        wjjs = [w for w in Collection(event,"GenWjj","nGenWjj")]
        ts = []
        if self.debug: print "Event %d:%d:%d " % (event.run,event.lumi,event.evt)
        for it,t in enumerate(gens):
          if abs(t.pdgId) != 6: continue
          iGenWjj = -1; ib = -1
          for iW, W in enumerate(wjjs):
            if W.iMom == it: 
                iGenWjj = iW
                break
          if iGenWjj == -1: continue
          for i,b in enumerate(gens):
            if abs(b.pdgId) == 5 and b.motherIndex == it: 
                ib = i; break
          if ib == -1: continue
          rj,drj = closest(b, jets)
          rl,drl = closest(b, leps)
          if drj < 0.4 and drj <= drl:
            rhad = 1; irb = jets.index(rj)
          elif drl < 0.4 and drl < drj:
            rhad = 0; irb = leps.index(rl)
          else:
            irb = -1; rhad = -1
          ret['nGenTopWjj'] += 1;
          ret['GenTopWjj_iGenWjj'].append(iGenWjj)
          ret['GenTopWjj_ib'].append(ib)
          ret['GenTopWjj_itop'].append(it)
          ret['GenTopWjj_rhad'].append(rhad)
          ret['GenTopWjj_irb'].append(irb)
          if self.debug:
            print "GenTop %d: iWjj %d (ij %d %d, rj %d %d), ib %d, had %d, rb %d, rjs %d" % (it,iGenWjj,W.ij1,W.ij2,W.irj1,W.irj2,ib,rhad,irb,(W.irj1 != -1) + (W.irj2 != -1))
        return ret

class WjjReco:
    def __init__(self,debug=False,extraVars=True,presel = lambda j1,j2,jjp4 : 30 < jjp4.M() and jjp4.M() < 130 and max(j1.btagCSV,j2.btagCSV) < 0.99 and min(j1.btagCSV,j2.btagCSV)<0.89): 
                                    #,presel=lambda j1,j2,jjp4 : True):
        self.debug = debug
        self.extravars = extraVars
        self.presel = presel
        self.branches = [ ("nWjj",'I'),
                 ("Wjj_pt","F",50,"nWjj"),
                 ("Wjj_eta","F",50,"nWjj"),
                 ("Wjj_phi","F",50,"nWjj"),
                 ("Wjj_mass","F",50,"nWjj"),
                 ("Wjj_ij1","I",50,"nWjj"),
                 ("Wjj_ij2","I",50,"nWjj"),
                 ("Wjj_mciW","I",50,"nWjj"),   
                 ("Wjj_mcij1","I",50,"nWjj"),
                 ("Wjj_mcij2","I",50,"nWjj") ]
        if self.extravars:
            self.branches += [
                 ("Wjj_j1pt","F",50,"nWjj"),
                 ("Wjj_j2pt","F",50,"nWjj"),
                 ("Wjj_dEta","F",50,"nWjj"),
                 ("Wjj_rap","F",50,"nWjj"),
                 ("Wjj_j1btagCSV","F",50,"nWjj"),
                 ("Wjj_j2btagCSV","F",50,"nWjj"),
                 ("Wjj_j1qgl","F",50,"nWjj"),
                 ("Wjj_j2qgl","F",50,"nWjj") ]
    def listBranches(self):
        return self.branches
    def __call__(self,event):
        # make python lists as Collection does not support indexing in slices
        jets = [j for j in Collection(event,"Jet","nJet")]
        # prepare output
        ret = { self.branches[0][0]: 0 }
        for b in self.branches[1:]: ret[b[0]] = []
        # do combinatorics
        for ij1,j1 in enumerate(jets):
            for ij2,j2 in enumerate(jets):
                if ij2 <= ij1: continue
                jjp4 = j1.p4() + j2.p4()
                if not self.presel(j1,j2,jjp4): continue
                ret["nWjj"] += 1;
                ret["Wjj_pt"].append(jjp4.Pt());
                ret["Wjj_eta"].append(jjp4.Eta());
                ret["Wjj_phi"].append(jjp4.Phi());
                ret["Wjj_mass"].append(jjp4.M());
                ret["Wjj_ij1"].append(ij1);
                ret["Wjj_ij2"].append(ij2);
                if self.extravars:
                    ret["Wjj_j1pt"].append(j1.pt)
                    ret["Wjj_j1btagCSV"].append(j1.btagCSV)
                    ret["Wjj_j1qgl"].append(j1.qgl)
                    ret["Wjj_j2pt"].append(j2.pt)
                    ret["Wjj_j2btagCSV"].append(j2.btagCSV)
                    ret["Wjj_j2qgl"].append(j2.qgl)
                    ret["Wjj_dEta"].append(abs(j1.eta-j2.eta))
                    ret["Wjj_rap"].append(jjp4.Rapidity())
                mcMatch = False
                if event.run == 1:
                    if self.debug: print "pair (%d,%d), M = %.2f" % (ij1,ij2,jjp4.M())
                    for i in xrange(event.nGenWjj):
                        gen_irj1 = event.GenWjj_irj1[i]
                        gen_irj2 = event.GenWjj_irj2[i]
                        if self.debug: print " ---> gen[%d] = (%d,%d)" % (i,gen_irj1,gen_irj2)
                        if (gen_irj1,gen_irj2) in [ (ij1,ij2), (ij2,ij1) ]:
                            if self.debug: print "        ---> yes "  
                            ret["Wjj_mciW"].append( event.GenWjj_iW[i] )
                            ret["Wjj_mcij1"].append( event.GenWjj_ij1[i] )
                            ret["Wjj_mcij2"].append( event.GenWjj_ij2[i] )
                            mcMatch = True;
                            break
                if not mcMatch:
                    ret["Wjj_mciW"].append( -1 )
                    ret["Wjj_mcij1"].append( -1 )
                    ret["Wjj_mcij2"].append( -1 )
        return ret

class WjjBestCand:
    def __init__(self, sortby = "mass", maxcand = 3, label = "", mva = None, debug=False):
        self.sortby = sortby
        self.maxcand = maxcand
        self.label = label
        self.debug = debug
        self.copyfloats = [ "pt", "eta", "phi", "mass" ]
        self.copyints = [ "ij1", "ij2", "mciW", "mcij1", "mcij2" ]
        self.branches = [ ("nBestWjj"+self.label,'I'), ("BestWjj"+self.label+"_index","I",self.maxcand,"nBestWjj"+self.label) ]
        self.mva = mva
        for v in self.copyfloats: 
            self.branches.append( ("BestWjj"+self.label+"_"+v,"F",self.maxcand,"nBestWjj"+self.label) )
        for v in self.copyints: 
            self.branches.append( ("BestWjj"+self.label+"_"+v,"I",self.maxcand,"nBestWjj"+self.label) )
        if self.mva: 
            self.branches.append( ("BestWjj"+self.label+"_"+"mva","F",self.maxcand,"nBestWjj"+self.label) )
    def listBranches(self):
        return self.branches
    def __call__(self,event):
        # make python lists as Collection does not support indexing in slices
        wjj = [ w for w in Collection(event,"Wjj","nWjj") ]
        for (i,w) in enumerate(wjj):
            w.index = i
        used = {}
        retlist = []
        if self.mva:
            for w in wjj: w.mva = self.mva(w)
        if self.sortby == "mass":
            wjj.sort(key = lambda w : abs(w.mass - 80.4))
        else:
            wjj.sort(key = lambda w : w.mva, reverse = True)
        if self.debug and sum(w.mciW != -1 for w in wjj) : 
            print "\n[%s]: Event %d:%d:%d with %d jets, %d Wjj candidates, of which %d mc matched " % (self.label, event.run,event.lumi,event.evt, event.nJet25, len(wjj), sum(w.mciW != -1 for w in wjj))
            used_printout = {}
            for w in wjj:
                locked = w.ij1 in used_printout or w.ij2 in used_printout
                print "[%s]\tcandidate %2d mass %6.1f   j1 %2d j2 %2d  score %+8.3f  mc match %-5s locked %-5s" % (
                             self.label, w.index, w.mass, w.ij1, w.ij2, w.mva if self.mva else abs(w.mass-80.4), w.mciW != -1, locked  )
                if not locked:
                    used_printout[w.ij1] = True
                    used_printout[w.ij2] = True
        for w in wjj:
            if w.ij1 in used or w.ij2 in used: continue
            retlist.append(w)
            used[w.ij1] = True
            used[w.ij2] = True
            if len(retlist) == self.maxcand: break
        # prepare output
        ret = { "nBestWjj"+self.label: len(retlist) }
        for v in self.copyfloats + self.copyints + ["index"]:
            ret[ "BestWjj"+self.label+"_"+v ] = [ getattr(w,v) for w in retlist ]       
        if self.mva: 
            ret[ "BestWjj"+self.label+"_mva" ] = [ w.mva for w in retlist ]       
        return ret
def WjjBestCandLikelihood(training,xml):
    from CMGTools.TTHAnalysis.tools.mvaTool import MVAVar, MVATool
    if training == "bestWjj":
        return MVATool("Likelihood",xml,vars = [
            MVAVar("Wjj_mass", func = lambda w : w.mass),
            MVAVar("Wjj_pt", func = lambda w : w.pt),
            MVAVar("Wjj_j2pt", func = lambda w : w.j2pt),
            MVAVar("Wjj_maxbtagCSV := max(0,max(Wjj_j1btagCSV,Wjj_j2btagCSV))", func = lambda w : max(0,max(w.j1btagCSV,w.j2btagCSV))),
            MVAVar("Wjj_minqgl := min(Wjj_j1qgl,Wjj_j2qgl)", func = lambda w : max(0,max(w.j1qgl,w.j2qgl)))
        ])

class TopWjjReco:
    def __init__(self,wjj="BestWjj", label = "", isBest=True, debug=False, extraVars=True, presel = lambda w,b,p4,had : True): 
        self.wjj = wjj
        self.label = label
        self.debug = debug
        self.isBest = isBest
        self.extravars = extraVars
        self.presel = presel
        self.branches = [ ("nTopWjj"+self.label,'I'),
                 ("TopWjj"+self.label+"_pt","F",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_eta","F",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_phi","F",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_mass","F",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_iW", "I",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_ib", "I",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_had","I",50,"nTopWjj"+self.label),
                 ("TopWjj"+self.label+"_mcitop", "I",50,"nTopWjj"+self.label) ]
    def listBranches(self):
        return self.branches
    def __call__(self,event):
        # make python lists as Collection does not support indexing in slices
        jets = [j for j in Collection(event,"Jet","nJet")]
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        wjjs = [j for j in Collection(event,self.wjj,"n"+self.wjj)]
        allwjjs = [j for j in Collection(event,"Wjj","nWjj")]
        if event.run == 1:
            genw = [j for j in Collection(event,"GenWjj",   "nGenWjj")   ]
            gent = [j for j in Collection(event,"GenTopWjj","nGenTopWjj")]
        # prepare output
        ret = { self.branches[0][0]: 0 }
        for b in self.branches[1:]: ret[b[0]] = []
        # do combinatorics
        if self.debug: print "Event %d:%d:%d " % (event.run,event.lumi,event.evt)
        for iw,w in enumerate(wjjs):
            thew = allwjjs[w.index] if self.isBest else w
            items =  [ (ib,b,1) for (ib,b) in enumerate(jets) ]
            items += [ (il,l,0) for (il,l) in enumerate(leps) ]
            for ib,b,had in items:
                if had == 1 and (ib == thew.ij1 or ib == thew.ij2): continue
                bp4 = b.p4() if had else b.p4() * (1.0/b.jetPtRatio)
                p4 = thew.p4() + bp4
                if not self.presel(thew,b,p4,had): continue
                ret["nTopWjj"+self.label] += 1;
                ret["TopWjj"+self.label+"_pt"].append(p4.Pt());
                ret["TopWjj"+self.label+"_eta"].append(p4.Eta());
                ret["TopWjj"+self.label+"_phi"].append(p4.Phi());
                ret["TopWjj"+self.label+"_mass"].append(p4.M());
                ret["TopWjj"+self.label+"_iW"].append(w.index if self.isBest else iw);
                ret["TopWjj"+self.label+"_ib"].append(ib);
                ret["TopWjj"+self.label+"_had"].append(had)
                mcitop = -1
                if event.run == 1:
                    for igt,gt in enumerate(gent):
                        if gt.rhad != had or gt.irb != ib: continue
                        gw = genw[gt.iGenWjj]
                        if thew.mciW != gw.iW: continue
                        mcitop = gt.itop
                ret["TopWjj"+self.label+"_mcitop"].append(mcitop)
                if self.debug:
                   print "RecTop: iWjj %2d (ij %d %d mc % 3d), ib %d, had %d, mass %7.1f, mc %s" % (w.index if self.isBest else iw, thew.ij1, thew.ij2, thew.mciW, ib, had, p4.M(), mcitop != -1)
        return ret


class TopReco_MC:
    def __init__(self):
        self.branches = [ ]
    def listBranches(self):
        return [ "mc_"+x for x in self.branches ]
    def __call__(self,event):
        ## prepare output container
        ret  = dict([(name,0.0) for name in self.branches])
        ## return if not on MC
        ret0 = dict([("mc_"+name,0.0) for name in self.branches])
        if event.run != 1:   return ret0
        #if event.nJet < 3: return ret0
        # make python lists as Collection does not support indexing in slices
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets = [j for j in Collection(event,"Jet","nJet")]
        gens = [j for j in Collection(event,"GenPart","nGenPart")]
        bjets = [ j for j in jets if j.btagCSV > 0.814 ]
        if len(bjets) == 0: bjets.append(jets[0])
        (met, metphi)  = event.met_pt, event.met_phi
        metp4 = ROOT.TLorentzVector()
        metp4.SetPtEtaPhiM(met,0,metphi,0)
        njet = len(jets); nb = len(bjets); nlep = len(leps)
        print "\nEvent %d:%d:%d " % (event.run,event.lumi,event.evt)
        for i,l in enumerate(leps):
            print "lepton %2d pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d   mcMatchId %+3d mcMatchAny %+3d " % ( 
                i, l.pt, l.eta, l.phi, l.pdgId, l.mcMatchId, l.mcMatchAny)
        for i,j in enumerate(jets):
            print "jet    %2d pT %7.2f eta %+5.2f phi %+5.2f btag %.3f mass %6.1f (pruned %6.1f)  mcMatchId %+3d mcFlavour  %+3d " % (
                i, j.pt, j.eta, j.phi, max(j.btagCSV,0.), j.mass, j.prunedMass, j.mcMatchId, j.mcFlavour)
        #for ig,g in enumerate(gens):
        #    print "gp %3d pT %7.2f eta %+5.2f phi %+5.2f pdgId %+3d    momIdx% 3d momId %+3d sourceId %+3d " % ( 
        #        ig, g.pt, g.eta, g.phi, g.pdgId, g.motherIndex, (g.motherId if g.motherId != -9999 else -99), g.sourceId)
        print ""
        tops = []
        for itop, top in enumerate(gens):
            if abs(top.pdgId) == 6:
                W = None; j1, j2, ij1, ij2 = None, None, -1, -1; l,v,il,iv = (None,None,None,None)
                for iW,cW in enumerate(gens):
                    if cW.motherIndex == itop and abs(cW.pdgId)==24:
                        W = cW; break
                if W:
                    for il,cl in enumerate(gens):
                        if cl.motherIndex == iW and abs(cl.pdgId)>6 and cl.charge != 0:
                            l = cl; break
                    for iv,cv in enumerate(gens):
                        if cv.motherIndex == iW and abs(cv.pdgId)>6 and cv.charge == 0:
                            v = cv; break
                    for ij1,cj1 in enumerate(gens):
                        if cj1.motherIndex == iW and abs(cj1.pdgId)<6:
                            for ij2, cj2 in enumerate(gens):
                                if ij2 <= ij1: continue
                                if cj2.motherIndex == iW and abs(cj2.pdgId)<6:
                                    j1 = cj1; 
                                    j2 = cj2; 
                                    if (j1.pt < j2.pt): 
                                        (j1,j2,ij1,ij2) = (j2,j1,ij2,ij1)
                                    break
                            break
                b = None; ib = None
                for ib,cb in enumerate(gens):
                    if cb.motherIndex == itop and abs(cb.pdgId)==5:
                        b = cb; break
                top.W = W; top.b = b; top.j1 = j1; top.j2 = j2; top.l = l; top.v = v
                top.idx = { 'top':itop, 'W':iW, 'j1':ij1, 'j2':ij2, 'l':il, 'v':iv, 'b':ib }
                tops.append(top) 
                if W:
                    if l:
                        rl,dr = closest(l, leps, lambda gen,rec: abs(gen.pdgId) == abs(rec.pdgId))
                        if dr < 1.0: l.rec = rl
                        else:        l.rec = None
                        rj,drj = closest(l, jets)
                        if drj < 0.5 and drj < dr: l.jrec = rj
                        else:                      l.jrec = None
                    if v:
                        v.rec = metp4
                    for (ij,j) in (ij1,j1), (ij2,j2):
                        if not j: continue
                        rj,dr = closest(j, jets)
                        if dr < 0.5: j.rec = rj
                        else:        j.rec = None
                    if j1 and j2 and j1.rec and j2.rec:
                        W.jj = j1.rec.p4() + j2.rec.p4()
                    else:
                        W.jj = None
                    W.mtlv = mt(l.rec, v.rec) if l and v and l.rec else 0
                    W.lv = l.rec.p4() + v.rec if l and v and l.rec else 0
                if b:
                    rj,dr = closest(b, jets)
                    if dr < 0.5: b.rec = rj
                    else:        b.rec = None
                    rl,drl = closest(b, leps)
                    if drl < 0.5 and drl < dr: b.lrec = rl
                    else:                      b.lrec = None
                top.jjb   = (W.jj +  b.rec.p4()) if W.jj and b and b.rec  else None
                top.jjlb  = (W.jj + jp4(b.lrec)) if W.jj and b and b.lrec else None
                top.sjjb  = (scalWjj(W.jj) +  b.rec.p4()) if W.jj and b and b.rec  else None
                top.sjjlb = (scalWjj(W.jj) + jp4(b.lrec)) if W.jj and b and b.lrec else None
                top.lvb   = toplvb(b.rec.p4(),l.rec,v.rec) if l and v and l.rec and b and b.rec  else None
                top.mtlvb  = mt(b.rec,  l.rec, v.rec) if l and v and l.rec and b and b.rec  else 0
                top.mtlvlb = mt(b.lrec, l.rec, v.rec) if l and v and l.rec and b and b.lrec else 0
        for top in tops:
                W = top.W; b = top.b; j1 = top.j1; j2 = top.j2; l = top.l; v = top.v
                itop = top.idx['top']; iW = top.idx['W']; ij1 = top.idx['j1']; ij2 = top.idx['j2']; il = top.idx['l']; iv = top.idx['v']; ib = top.idx['b']
                print "top %2d       pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( itop, top.pt, top.eta, top.phi, top.pdgId)
                if top.jjb:
                    print " |      jjb  pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                                top.jjb.Pt(), top.jjb.Eta(), top.jjb.Phi(), top.jjb.M(), deltaR(top.jjb.Eta(),top.jjb.Phi(),top.eta,top.phi), top.jjb.Pt()/top.pt )
                    print " |     sjjb  pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                                top.sjjb.Pt(), top.sjjb.Eta(), top.sjjb.Phi(), top.sjjb.M(), deltaR(top.sjjb.Eta(),top.sjjb.Phi(),top.eta,top.phi), top.sjjb.Pt()/top.pt )
                if top.lvb:
                    print " |      lvb  pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f mt %5.1f   dr %5.3f ptR %4.2f  " % ( 
                                top.lvb.Pt(), top.lvb.Eta(), top.lvb.Phi(), top.lvb.M(), top.mtlvb, deltaR(top.lvb.Eta(),top.lvb.Phi(),top.eta,top.phi), top.lvb.Pt()/top.pt )
                if top.jjlb:
                    print " |      jjlb pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                                top.jjlb.Pt(), top.jjlb.Eta(), top.jjlb.Phi(), top.jjlb.M(), deltaR(top.jjlb.Eta(),top.jjlb.Phi(),top.eta,top.phi), top.jjlb.Pt()/top.pt )
                    print " |     sjjlb pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                                top.sjjlb.Pt(), top.sjjlb.Eta(), top.sjjlb.Phi(), top.sjjlb.M(), deltaR(top.sjjlb.Eta(),top.sjjlb.Phi(),top.eta,top.phi), top.sjjlb.Pt()/top.pt )
                print " | "
                if W:
                    print " |- W  %2d    pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( iW, W.pt, W.eta, W.phi, W.pdgId)
                    if j1 and j2 and j1.rec and j2.rec:
                        print " |  |     jj pT %7.2f eta %+5.2f phi %+5.2f mass %5.1f            dr %5.3f ptR %4.2f  " % ( 
                                W.jj.Pt(), W.jj.Eta(), W.jj.Phi(), W.jj.M(), deltaR(W.jj.Eta(),W.jj.Phi(),W.eta,W.phi), W.jj.Pt()/W.pt )
                    if W.lv:
                        print " |  |     lv pT %7.2f           phi %+5.2f mT   %5.1f            dp %5.3f ptR %4.2f  " % ( 
                                W.lv.Pt(),             W.lv.Phi(), W.mtlv, abs(deltaPhi(W.lv.Phi(),W.phi)), W.lv.Pt()/W.pt )
                    if l:
                        print " |  +- l %2d  pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( il, l.pt, l.eta, l.phi, l.pdgId)
                        if l.rec:
                            print " |  |    lep pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d btag %.3f dr %5.3f ptR %4.2f ptRj %4.2f " % ( 
                                l.rec.pt, l.rec.eta, l.rec.phi, l.rec.pdgId, max(l.rec.jetBTagCSV,0), deltaR(l,l.rec), l.rec.pt/l.pt, l.rec.pt/l.rec.jetPtRatio/l.pt )
                        if l.jrec:
                            print " |  |    jet pT %7.2f eta %+5.2f phi %+5.2f            btag %.3f dr %5.3f ptR %4.2f" % ( 
                                l.jrec.pt, l.jrec.eta, l.jrec.phi, max(0,l.jrec.btagCSV), deltaR(l,l.jrec), l.jrec.pt/l.pt )
                    print " |  |  " 
                    if v:
                        print " |  +- v %2d  pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( iv, v.pt, v.eta, v.phi, v.pdgId)
                        print " |       met pT %7.2f           phi %+5.2f                       dp %5.3f ptR %4.2f" % ( 
                                        v.rec.Pt(), v.rec.Phi(), abs(deltaPhi(v.phi,v.rec.Phi())), v.rec.Pt()/v.pt )
                    for (ij,j) in (ij1,j1), (ij2,j2):
                        if not j: continue
                        if ij == ij2: print " |  |  " 
                        print " |  +- j %2d  pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( ij, j.pt, j.eta, j.phi, j.pdgId)
                        if j.rec:
                            print " |  %s    jet pT %7.2f eta %+5.2f phi %+5.2f            btag %.3f dr %5.3f ptR %4.2f" % ( " " if ij == ij2 else "|",
                                j.rec.pt, j.rec.eta, j.rec.phi, max(0,j.rec.btagCSV), deltaR(j,j.rec), j.rec.pt/j.pt )
                else:
                    print " |- W     not found!!"
                print " | "
                if b:
                    print " +- b %2d     pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d " % ( ib, b.pt, b.eta, b.phi, b.pdgId)
                    if b.rec:
                        print "         jet pT %7.2f eta %+5.2f phi %+5.2f            btag %.3f dr %5.3f ptR %4.2f" % ( 
                            b.rec.pt, b.rec.eta, b.rec.phi, max(0,b.rec.btagCSV), deltaR(b,b.rec), b.rec.pt/b.pt )
                    if b.lrec:
                        print "         lep pT %7.2f eta %+5.2f phi %+5.2f pdgId  %+3d btag %.3f dr %5.3f ptR %4.2f ptRj %4.2f " % ( 
                            b.lrec.pt, b.lrec.eta, b.lrec.phi, b.lrec.pdgId, max(b.lrec.jetBTagCSV,0), deltaR(b,b.lrec), b.lrec.pt/b.pt, b.lrec.pt/b.lrec.jetPtRatio/b.pt )
                else:
                    print " +- b      not found!!"

                print ""
        print "-" * 100
        print ""

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.truth  = TopReco_MC()
            self.mc  = GenWjj(debug=False)
            self.mct  = GenTopWjj(debug=True)
            self.rec = WjjReco(debug=False,extraVars=True)
            self.mvasort = WjjBestCandLikelihood("bestWjj","/afs/cern.ch/work/g/gpetrucc/CMSSW_7_4_7/src/CMGTools/TTHAnalysis/weights/bestWjj_Likelihood.weights.xml")
            self.sort = WjjBestCand(debug=True, mva = self.mvasort, sortby = "mva")
            self.msort = WjjBestCand(debug=False, label = "ByMass")
            self.top = TopWjjReco(debug=True, presel = lambda w,b,p4,had : had != 1 or b.btagCSV > 0.89)
            self.njnb = [    ((3,3),(1,1)), ((3,3),(2,2)), ((3,3),(3,999)),
                             ((4,4),(1,1)), ((4,4),(2,2)), ((4,4),(3,999)),
                             ((5,5),(1,1)), ((5,5),(2,2)), ((5,5),(3,999)),
                             ((6,999),(1,1)), ((6,999),(2,2)), ((6,999),(3,999)) ]
            self._acc = {} ; self._obs = {} ; self._sort = {} ; self._msort = {}
            for (nj,nb) in self.njnb:
                self._acc[(nj,nb)] = [0,0,0,0,0,0]
                self._obs[(nj,nb)] = [0,0,0]
                self._sort[(nj,nb)] = [[0,0,0],[0,0,0],[0,0,0]]
                self._msort[(nj,nb)] = [[0,0,0],[0,0,0],[0,0,0]]
        def analyze(self,ev):
            mybin = None
            ret = self.mc(ev)
            for abin in self.njnb:
                ((njmin,njmax),(nbmin,nbmax)) = abin
                if njmin <= ev.nJet and ev.nJet <= njmax and nbmin <= ev.nBJetMedium25 and ev.nBJetMedium25 <= nbmax:
                    mybin = abin
                    break
            if not mybin: return False
            self._acc[mybin][0] += 1
            ret = self.mc(ev)
            self._acc[mybin][1] += ret['nGenWjj']
            found = []
            for i in xrange(ret['nGenWjj']):
                if ret['GenWjj_acc'][i]: 
                    self._acc[mybin][2] += 1
                    if ret['GenWjj_irj1'][i] != -1 and ret['GenWjj_irj2'][i] != -1: 
                        self._acc[mybin][3] += 1
                        self._acc[mybin][4] += ( ret['GenWjj_irj1'][i] != ret['GenWjj_irj2'][i] ) 
                elif ret['GenWjj_irj1'][i] != -1 and ret['GenWjj_irj2'][i] != -1:
                    if ret['GenWjj_irj1'][i] != -1 and ret['GenWjj_irj2'][i] != -1: 
                        self._acc[mybin][5] += 1
            for k,v in ret.iteritems(): setattr(ev,k,v)
            rett = self.mct(ev)
            for k,v in rett.iteritems(): setattr(ev,k,v)
            if rett["nGenTopWjj"]: self.truth(ev)
            ret2 = self.rec(ev)
            self._obs[mybin][0] += 1
            #if self._obs[0] < 20: print ret2
            self._obs[mybin][1] += ret2['nWjj']
            self._obs[mybin][2] += len([ i for i in xrange(ret2['nWjj']) if ret2['Wjj_mciW'][i] != -1 ])
            for k,v in ret2.iteritems():
                setattr(ev,k,v)
            ret3 = self.sort(ev)
            ret3m = self.msort(ev)
            for k,v in ret3.iteritems(): setattr(ev,k,v)
            for k,v in ret3m.iteritems(): setattr(ev,k,v)
            ret4 = self.top(ev)
            nacc = 0; nrecmatch = 0
            for i in xrange(ret['nGenWjj']):
                if not ret['GenWjj_acc'][i]: continue
                nacc += 1
            nrecmatch = sum(x != -1 for x in ret2["Wjj_mciW"])
            for i in xrange(3):
                if nacc > i:
                    self._sort[mybin][i][0] += 1
                    self._msort[mybin][i][0] += 1
                if nrecmatch > i:
                    self._sort[mybin][i][1] += 1
                    self._msort[mybin][i][1] += 1
                    if ret3["nBestWjj"] > i and ret3["BestWjj_mciW"][i] != -1:
                        self._sort[mybin][i][2] += 1
                    if ret3m["nBestWjjByMass"] > i and ret3m["BestWjjByMass_mciW"][i] != -1:
                        self._msort[mybin][i][2] += 1
    test = Tester("tester")              
    el = EventLoop([ test ])
    el.loop([tree], maxEvents = 100000 if len(argv) <= 2 else int(argv[2]))
    if False:
      for abin in test.njnb:
        ((njmin,njmax),(nbmin,nbmax)) = abin
        print "bin for %d <= nJet <= %d, %d <= nBJet <= %d" % (njmin,njmax,nbmin,nbmax)
        print "\t",test._acc[abin]
        print "\t",test._obs[abin]
        print "\t",test._sort[abin]
        print "\t",test._msort[abin]
        print
