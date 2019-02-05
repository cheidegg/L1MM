import types
import time
import ROOT
import pandas

from functions           import *
from functionsTreeReader import *


class Event:
	def __init__(self, master, tree, entry):
		self.master = master
		self.tree   = tree
		self.entry  = entry
		self.objs   = {}
		self.load()
	def __getattr__(self, key):
		if key in self.objs: return self.objs[key]
		#if hasattr(self.tree, key): return getattr(self.tree, key)
		return readBranch(self.tree, key)
	def __getitem__(self, attr):
		return self.__getattr__(attr)
	def getVector(self, key, index):
		return readBranch(self.tree, key, index)
	def load(self):
		if self.tree.entry == self.entry: return
		if self.tree.entry == self.entry-1: self.tree._ttreereader.Next()
		else: self.tree._ttreereader.SetEntry(self.entry)
		self.tree.entry = self.entry



class Object:
	def __init__(self, master, event, objDef, index):
		self.master = master
		self.event  = event
		self.objDef = objDef
		self.index  = index
		self.values = {}
	def __eq__(self, other):
		if other.__class__.__name__!="Object": return False
		if self is other: return True
		if self.Et==other.Et and self.Eta==other.Eta and self.Phi==other.Phi: return True
		return False
	def __ne__(self, other):
		if other.__class__.__name__!="Object": return True
		if self is other: return False
		if self.Et==other.Et and self.Eta==other.Eta and self.Phi==other.Phi: return False
		return True
	def __getattr__(self, key):
		if key in self.values.keys(): return self.values[key]
		value = self.event.getVector(self.objDef.branches[key], self.index)
		self.values[key] = value
		return value
	def __getitem__(self, attr):
		return self.__getattr__(attr)
	def p4(self, pt = 0):
		ret = ROOT.TLorentzVector()
		if pt > 0: ret.SetPtEtaPhiM(pt     ,self.eta,self.phi,self.mass)
		else     : ret.SetPtEtaPhiM(self.pt,self.eta,self.phi,self.mass)
		return ret
	def __repr__(self):
		return ("<%s[%s]>" % (self.objDef.name, self.index) if self.index != None else ("<%s>" % self.objDef.name))
	def __str__(self):
		return self.__repr__()



class Objectlist:
	def __init__(self, master, event, objDef):
		self.master = master
		self.event  = event
		self.objDef = objDef
		self.name   = objDef.name
		if objDef.isFlat:
			self.length = 1
		else:
			self.length = int(getattr(event, objDef.lengthbranch))
			self.length = self.length if self.length else 0
		self.objs   = {}
	def __getitem__(self, index):
		if type(index) == int and index in self.objs: return self.objs[index]
		if index >= self.length: self.master.vb.error("Invalid index %r (length is %r) for object %s" % (index, self.length, self.objDef.name))
		obj = Object(self.master, self.event, self.objDef, index)
		self.objs[index] = obj
		return obj
	def __len__(self):
		return self.length


class Sample(object):
	## mother class for Samples
	def __init__(self, master, sampleDef):
		self.master  = master
		self.name    = sampleDef.name
		self.opts    = sampleDef.options
		self.opened  = False
		self.characterize()
	def characterize(self):
		## find and store the virtual id of the sample
		self.sampleId = self.master.getSampleId(self.opts["path"])
	def close(self):
		## close the ntuple
		if not self.opened: return
		self.tfile.Close()
		self.opened = False	
	def getEntries(self):
		## return total entries in the ntuple
		if not hasattr(self, "entries"): self.setEntries()
		return self.entries
	def init(self):
		## initializing the tree reader thing
		self.tree.entry = -1
		self.tree._ttreereader = ROOT.TTreeReader(self.tree)
		self.tree._ttreereader.SetEntry(0)
		self.tree._ttrvs = {}
		self.tree._ttras = {}
		self.tree._leafTypes = {}
		self.tree._ttreereaderversion = 1
		self.tree.arrayReader = types.MethodType(getArrayReader, self.tree)
		self.tree.valueReader = types.MethodType(getValueReader, self.tree)
		self.tree.readBranch  = types.MethodType(readBranch    , self.tree)
	def load(self):
		self.open()
		self.init()
	def open(self):
		## open the ntuple for the sample
		if self.opened: return
		self.tfile  = ROOT.TFile.Open(self.opts["path"],"read")
		self.tree   = getRootObj(self.tfile, self.opts["tree"])
		self.opened = True
	def setEntries(self):
		if not self.opened: self.load()
		self.entries = self.tree.GetEntries()

class SampleDefault(Sample):
	## default sample class using TEventLists
	def __init__(self, master, sampleDef):
		super(SampleDefault, self).__init__(master, sampleDef)
		self.lists = {} # ROOT TEventList
		self.retrieveLists()
	def analyze(self, triggers):
		## building the TEventLists for every trigger
		if not self.opened: self.load()
		self.setEntries()
		toRun = []
		for trigger in triggers:
			if trigger.triggerId in self.lists  .keys() and not self.master.getOpt("force"): continue
			self.lists  [trigger.triggerId] = ROOT.TEventList()
			toRun.append(trigger.triggerId)
		if len(toRun)==0: return
		t0       = time.clock()
		tlast    = t0
		for ievt in xrange(self.getEntries()):
			event       = Event(self.master, self.tree, ievt)
			trigObjects = {}
			for objDef in self.master.menu.objects:
				trigObjects[objDef.name] = Objectlist(self.master, event, objDef)
			for it,trigger in enumerate([x for x in triggers if x.triggerId in toRun]):
				if trigger.apply(event, trigObjects): 
					self.lists[trigger.triggerId].Enter(ievt)
			if ievt>0 and ievt%10000 == 0:
				t1 = time.clock()
				self.master.vb.talk("Processed %8d/%8d entries of this tree (elapsed time %7.1fs, curr speed %8.3f kHz, avg speed %8.3f kHz)" % (ievt, self.entries, t1-t0, (10.000)/(max(t1-tlast,1e-9)),ievt/1000./(max(t1-t0,1e-9))))
				tlast = t1
		for newId in toRun:
			self.master.addList(newId, self.sampleId, self.lists[newId])
	def apply(self, trigger):
		## run a single trigger
		if not trigger.triggerId in self.lists.keys(): 
			self.analyze([trigger])
		return self.lists[trigger.triggerId].GetN()
	def intersect(self, triggers):
		## building the list for events passing both triggers together
		toRun = []
		for trigger in triggers:
			if trigger.triggerId in self.lists: continue
			toRun.append(trigger)
		if len(toRun)>0: self.analyze(toRun)
		newId = self.master.getTrigIdFromList(triggers, False)
		lists = [self.lists[trigger.triggerId] for trigger in triggers]
		intersected = lists[0].Clone(newId)
		for thelist in lists[1:]:
			intersected.Intersect(thelist)
		self.lists[newId] = intersected
		self.master.addList(newId, self.sampleId, intersected)
		return intersected.GetN()
	def merge(self, triggers):
		## building the list for events passing any of the triggers
		toRun = []
		for trigger in triggers:
			if trigger.triggerId in self.lists: continue
			toRun.append(trigger)
		if len(toRun)>0: self.analyze(toRun)
		newId = self.master.getTrigIdFromList(triggers, True)
		lists = [self.lists[trigger.triggerId] for trigger in triggers]
		merged = lists[0].Clone(newId)
		for thelist in lists[1:]:
			merged.Add(thelist)
		self.lists[newId] = merged
		self.master.addList(newId, self.sampleId, merged)
		return merged.GetN()
	def retrieveLists(self):
		## find all virtual trigger ids that contain this sample path
		## then compile a list of teventlists that apply to this sample
		self.lists = self.master.getListsToRetrieve(self.sampleId)



class SampleBuffers(Sample):
	## Sample class for use of Pandas DataFrame buffers
	def __init__(self, master, sampleDef):
		super(SampleBuffers, self).__init__(master, sampleDef)
		self.buffers = {} # Pandas DataFrame without selection
		self.sels    = {} # Pandas DataFrame with selection applied
		self.retrieveBuffers()
	def analyze(self, triggers):
		## building the pandas DataFrames for every trigger (not yet working)
		## the logic is the following:
		## * take the objects of a given trigger across all legs
		## * remove the threshold cut from all trigger legs
		## * build all pairs of all objects in the event 
		## * check whether they pass the trigger legs without threshold cut
		if not self.opened: self.load()
		self.setEntries()
		toRun     = []
		container = {}
		for trigger in triggers:
			if trigger.triggerIdeff in self.buffers.keys(): continue
			container[trigger.triggerIdeff] = []
			trigger.setThresholdCuts() ### ==> the thing is, i need to collect all threshold c
			toRun.append(trigger)
		if len(toRun)==0: return
		t0       = time.clock()
		tlast    = t0
		for ievt in xrange(self.getEntries()):
			event       = Event(self.master, self.tree, ievt)
			trigObjects = {}
			for objDef in self.master.menu.objects:
				trigObjects[objDef.name] = Objectlist(self.master, event, objDef)
			for it,trigger in enumerate(toRun):
				pairs = trigger.applyPairing(event, trigObjects)
				if len(pairs)==0: continue
				for ip, pair in enumerate(pairs):
					container[trigger.triggerIdeff].append([ievt, ip]+pair)
			if ievt>0 and ievt%10000 == 0:
				t1 = time.clock()
				self.master.vb.talk("Processed %8d/%8d entries of this tree (elapsed time %7.1fs, curr speed %8.3f kHz, avg speed %8.3f kHz)" % (ievt, self.entries, t1-t0, (10.000)/(max(t1-tlast,1e-9)),ievt/1000./(max(t1-t0,1e-9))))
				tlast = t1
		for it,trigger in enumerate(toRun):
			self.buffers[trigger.triggerIdeff] = pandas.DataFrame(container[trigger.triggerIdeff], columns=["ievt", "ipair"]+trigger.getFakedLegNames())
			self.master.addBuffer(trigger.triggerIdeff, self.sampleId, self.buffers[trigger.triggerIdeff])
	def apply(self, trigger):
		## evaluate the trigger based on the paired dataframe
		## note: we need to take the data frame that has been built with
		## the oldest parent, all children only have the thresholds varied
		## but the rest of the cuts remains the same, which is how the 
		## dataframe is built; this in fact is the key of doing this: we
		## take the combinations evaluated effortfully earlier and now 
		## only change the thresholds
		## the DF selection looks for events where at least one pair survives the cuts
		oldestParent = trigger.getOldestParent()
		if not oldestParent.triggerIdeff in self.buffers.keys(): 
		        self.analyze([oldestParent])
		df       = self.buffers[oldestParent.triggerIdeff] ## the data frame
		legnames = trigger.getFakedLegNames()
		#legops   = trigger.getFakedLegOps()
		legops   = [">" for i in range(len(legnames))] ## FIXME: should be derived from the cfg
		legcuts  = trigger.getFakedLegCuts() ### what are the required leg thresholds?
		sels     = ["(df."+legnames[il]+legops[il]+str(legcuts[il])+")" for il in range(len(legnames))]
		events   = eval("df.loc["+" & ".join(sels)+"].values.T.tolist()[0]")
		self.sels[trigger.triggerId] = list(set(events))
		return len(self.sels[trigger.triggerId])
	def createBuffer(self, triggers):
		## creating the pandas data frames, basically calling analyze
		if len(triggers)==0: return
		self.analyze(triggers)
	def intersect(self, triggers):
		## building the buffer for events passing both triggers together
		parents = [trigger.getOldestParent() for trigger in triggers]
		toRun = []
		for trigger in triggers:
			oldestParent = trigger.getOldestParent()
			if oldestParent.triggerIdeff in self.buffers.keys(): continue
			toRun.append(oldestParent)
		if len(toRun)>0: self.analyze(toRun)
		toRun = []
		for trigger in triggers:
			if trigger.triggerId in self.sels.keys(): continue
			toRun.append(trigger)
		if len(toRun)>0: 
			for trigger in toRun: self.apply(trigger)
		newId   = self.master.getTrigIdFromList(parents, False, True)
		buffers = [self.buffers[trigger.triggerIdeff] for trigger in parents ]
		sels    = [self.sels   [trigger.triggerId   ] for trigger in triggers]
		subsets = [subsetDF(buffers[idx], sels[idx]) for idx in range(len(buffers))]
		intersected, selected  = intersectDFs(buffers, subsets, ["ievt", "ipair"])
		self.buffers[newId] = intersected
		self.master.addBuffer(newId, self.sampleId, intersected)
		return len(list(set(selected["ievt"])))
	def merge(self, triggers):
		## building the list for events passing any of the triggers
		parents = [trigger.getOldestParent() for trigger in triggers]
		toRun = []
		for trigger in triggers:
			oldestParent = trigger.getOldestParent()
			if oldestParent.triggerIdeff in self.buffers.keys(): continue
			toRun.append(oldestParent)
		if len(toRun)>0: self.analyze(toRun)
		toRun = []
		for trigger in triggers:
			if trigger.triggerId in self.sels.keys(): continue
			toRun.append(trigger)
		if len(toRun)>0: 
			for trigger in toRun: self.apply(trigger)
		newId   = self.master.getTrigIdFromList(parents, True, True)
		buffers = [self.buffers[trigger.triggerIdeff] for trigger in parents ]
		sels    = [self.sels   [trigger.triggerId   ] for trigger in triggers]
		subsets = [subsetDF(buffers[idx], sels[idx]) for idx in range(len(buffers))]
		merged, selected  = mergeDFs(buffers, subsets, ["ievt", "ipair"])
		self.buffers[newId] = merged
		self.master.addBuffer(newId, self.sampleId, merged)
		return len(list(set(selected["ievt"])))
	def retrieveBuffers(self):
		## find all virtual trigger ids that contain this sample path
		## then compile a list of dataframes that apply to this sample
		self.buffers = self.master.getBuffersToRetrieve(self.sampleId)










