from importlib import import_module
class VIDWrapper(object):
    def __init__(self,factory,package):
        self._factory = factory
        self._package = package
        self._cache = {}
        self._defs  = {}
    def define(self,name,idname=None,cff=None,package=None,cutsToIgnore=None):
        if name in self._defs: raise RuntimeError, "%s already defined" % name
        if not idname:  idname = name
        if not package: package = self._package
        if not cff: raise RuntimeError, "cff must be provided in define() call"
        self._defs[name] = { 'idname':idname, 'cff':cff, 'package':package, 'cutsToIgnore':(cutsToIgnore[:] if cutsToIgnore else []) }
    def load(self,name,idname=None,cff=None,package=None,cutsToIgnore=None):
        if name not in self._cache:
            if name in self._defs:
                mydef = self._defs[name]
                if idname and (idname != mydef['idname']): raise RuntimeError, "ID name %s specified but not matching with definition for %s" % (idname,name)
                if cff and (cff != mydef['cff']): raise RuntimeError, "CFF %s specified but not matching with definition for %s" % (cff,name)
                if package and (package != mydef['package']): raise RuntimeError, "Package %s specified but not matching with definition for %s" % (package,name)
                if cutsToIgnore and (cutsToIgnore != mydef['cutsToIgnore']): raise RuntimeError, "Cuts to ignore %s specified but not matching with definition for %s" % (cutsToIgnore,name)
                idname, cff, package, cutsToIgnore = mydef['idname'], mydef['cff'], mydef['package'], mydef['cutsToIgnore']
                mod = import_module(cff,package)
            else:
                if not idname:  idname = name
                if not package: package = self._package
                if not cff:
                    if ":" in idname:
                        (fname,wp) = idname.split(":")
                        cff, idname = "."+fname+"_cff", fname+"_"+wp
                        mod = import_module(cff,package)
                    else:
                        try:
                            cff = "."+idname+"_cff"
                            mod = import_module(cff,package)
                        except ImportError:
                            (base,wp) = idname.rsplit("_",1)
                            cff, idname = "."+base+"_cff", base+"_"+wp
                            mod = import_module(cff,package)
                else:
                    mod = import_module(cff,package)
                    self._defs[name] = { 'idname':idname, 'cff':cff, 'package':package, 'cutsToIgnore':(cutsToIgnore[:] if cutsToIgnore else []) }
            pset = getattr(mod,idname)
            if cutsToIgnore:
                nfound = 0
                print name, cutsToIgnore
                for i,p in enumerate(pset.cutFlow):
                    print "\t"+p.cutName.value()
                    if p.cutName.value() in cutsToIgnore:
                        pset.cutFlow[i].isIgnored = True
                        nfound += 1
                if nfound != len(cutsToIgnore):
                    raise RuntimeError, "Error, at least one cut in %s was not found in %s" % (cutsToIgnore, pset)
            selector = self._factory(pset)
        self._cache[name] = selector
    def __getitem__(self,name,idname=None,cff=None,package=None,cutsToIgnore=None):
        if name not in self._cache:
            self.load(name,idname=idname,cff=cff,package=package,cutsToIgnore=cutsToIgnore)
        return self._cache[name]

from RecoEgamma.ElectronIdentification.VIDElectronSelector import VIDElectronSelector
from RecoEgamma.PhotonIdentification.VIDPhotonSelector import VIDPhotonSelector
from RecoMuon.MuonIdentification.VIDMuonSelector import VIDMuonSelector

ElectronVID = VIDWrapper(VIDElectronSelector, "RecoEgamma.ElectronIdentification.Identification")
PhotonVID = VIDWrapper(VIDPhotonSelector,"RecoEgamma.PhotonIdentification.Identification")
MuonVID = VIDWrapper(VIDMuonSelector, "RecoMuon.MuonIdentification.Identification")

ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_veto", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_veto", cff=".cutBasedElectronID_Spring15_25ns_V1_cff")
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_veto_noIso", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_veto", cff=".cutBasedElectronID_Spring15_25ns_V1_cff", cutsToIgnore=['GsfEleEffAreaPFIsoCut'])
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_loose", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_loose", cff=".cutBasedElectronID_Spring15_25ns_V1_cff")
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_loose_noIso", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_loose", cff=".cutBasedElectronID_Spring15_25ns_V1_cff", cutsToIgnore=['GsfEleEffAreaPFIsoCut'])
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_medium", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_medium", cff=".cutBasedElectronID_Spring15_25ns_V1_cff")
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_medium_noIso", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_medium", cff=".cutBasedElectronID_Spring15_25ns_V1_cff", cutsToIgnore=['GsfEleEffAreaPFIsoCut'])
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_tight", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_tight", cff=".cutBasedElectronID_Spring15_25ns_V1_cff")
ElectronVID.define(name="cutBasedElectronID_Spring15_25ns_V1_tight_noIso", idname="cutBasedElectronID_Spring15_25ns_V1_standalone_tight", cff=".cutBasedElectronID_Spring15_25ns_V1_cff", cutsToIgnore=['GsfEleEffAreaPFIsoCut'])

