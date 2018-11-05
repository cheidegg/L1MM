## this is the place to put hard-coded functions that can be
## used just like the ones defined on-the-fly in the cfg
## via the keyword "functions"; however, do not over-use
## this part, the on-the-fly definitions should be preferred
## whenever possible

import ROOT
import inspect

class HcFunctions:
	def __init__(self, master):
		self.master = master
	def _listAllFunctions(self):
		if hasattr(self, "fnames"): return self.fnames
		self.fnames = [x[0] for x in inspect.getmembers(self, predicate=inspect.ismethod)]
		return [x for x in self.fnames if x[0:1]!="_"]
	def _getFunctionEntity(self, functionName, *args, **kwargs):
		if not functionName in self._listAllFunctions(): return None
		return getattr(self, functionName)
	## --- put functions below ---
	def mass(self, objects):
		object1 = objects[0]; object2 = objects[1]
		o1 = ROOT.TLorentzVector()
		o1.SetPtEtaPhiM(object1.Et, object1.Eta, object1.Phi, object1.M if hasattr(object1, "M") else 0)
		o2 = ROOT.TLorentzVector()
		o2.SetPtEtaPhiM(object2.Et, object2.Eta, object2.Phi, object2.M if hasattr(object2, "M") else 0)
		return (o1+o2).M()


