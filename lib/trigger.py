from functions import *
from cfg       import *

import ROOT
import numpy
import math
import collections
import inspect

import functionsHardcoded

def _getArgs(trigLeg, selectedByOtherLegs):
	## builds the "args" collection, i.e. a list of objects that have been
	## selected by those previously evaluated trigger legs that are required
	## by this particular leg (if this leg does not make reference to objects
	## of the previous leg, the collection will be empty regardless of what is
	## in selectedByOtherLegs)
	args = [[] for i in range(len(trigLeg.functs.keys()))]
	for iu in range(len(args)):
		selected = [olist for io,olist in enumerate(selectedByOtherLegs) if io+1 in trigLeg.uses[iu]] ## first leg is called "1" not "0"
		nElms    = [len(theobjs) for theobjs in selected]
		nVars    = int(numpy.prod(nElms)) if len(nElms)>0 else 0
		if nVars == 0: continue
		for iVar in range(nVars):
			args[iu].append([selected[i][(iVar/numpy.prod(nElms[i+1:]))%numpy.prod(nElms[i])] for i in range(len(selected[0:-1]))] + [selected[-1][iVar%nElms[-1]]])
	return args



def _leg(functions, event, trigObjlists, trigLeg, selectedByOtherLegs, onlyThresholdCut=False, exceptThresholdCut=False):
	## this function evaluates a normal trigger leg (TriggerLeg) on an event that
	## contains a certain list of trigger objects (trigObjlists)
	## careful with overlapping legs: there can be multiple previously selected objects 
	## per leg => only need to find one of them to match; in fact, it would be better 
	## to have a loop over the different legs; however, this is time-consuming and only
	## really worth it if there are at least 3 legs of the same underlying object; for this
	## case, have a look at _multi and TriggerMulti; 
	## in all other cases, note: when running the second leg (that, let's assume so, depends 
	## on the first leg), the first leg is already run; it means, the possible overlapping
	## objects are already there in selectedByOtherLegs; thus we need to filter them out..

	result = []

	## collecting the trigger objects
	objs = trigObjlists[trigLeg.obj.name]

	## collecting the objects selected by previously evaluated trigger legs
	args = _getArgs(trigLeg, selectedByOtherLegs)

	## running the selection in the trigger leg per object
	## for this, pargs is extracted from the args list first
	selected = []
	for obj in [objs[io] for io in range(objs.length)]:
		pargs = []
		for iu in range(len(args)): 
			pargs.append([])
			for iv in range(len(args[iu])):
				pargs[iu].append([])
				pargs[iu][iv] = [x for x in args[iu][iv] if x!=obj]

		if not trigLeg.apply(functions, obj, pargs, onlyThresholdCut, exceptThresholdCut): 
			continue
		selected.append(obj)
	
	## removing any of the selected objects from the "selected" list in case
	## it is used in one of the previous legs already
	for obj in selected:
		used=False
		for io,others in enumerate(selectedByOtherLegs):
			if obj in others:
				if len(others)==1:   ## in this case, the other trigger leg 
					used=True; break ## would end up with len(objects)==0,
				                     ## which means the event is rejected
				                     ## so the event cannot be removed from
				                     ## that list nor can it be kept for this leg
				others.remove(obj)
		if not used: result.append(obj)
	return result



def _multiPermute(functions, functionCall, num, req, objlist): 
	## is called from within _multi (see below)
	## this function accounts for the fact that there can be n objects in 
	## the objlist collection (= combination built by the _multi selection)
	## among who a set of m < n objects enter the function call (string)
	## which implements an extra selection once the combination has been built

	def iterate(combination, maxlen):
		## generalized version of doing this here:
		## if combination[-1]<maxlen-1: combination[-1]+=1; return combination
		## if combination[-2]<maxlen-1: combination[-2]+=1; combination[-1]=1; return combination
		for j in reversed(range(len(combination))):
			combination[j] += 1
			if combination[j]==maxlen:
				combination[j] = combination[j-1]+2 ## the lower index will be increased 
				                                    ## by 1 in the next iteration, and we 
				                                    ## do not want two times the same
			else:
				return combination
		return combination

	## building the combinations of objects to enter the extra functions
	combinations = []
	last = [0 for i in range(num)]
	for ic in range(binomial(len(objlist), num)):
		## this loop maps the index ic to the object combinations
		## e.g. if there are 4 objects and we build pairs, we get
		## 6 combinations and the mapping looks like this:
		##   0 => 0 1 
		##   1 => 0 2
		##   2 => 0 3
		##   3 => 1 2
		##   4 => 1 3
 		##   5 => 2 3
		last    = iterate(last[:], len(objlist))
		combinations.append(last)

	## prepare the statement to be executed (the string functionCall is a 
	## statement, not only a function call, but it contains the function 
	## call)
	remastered = functionCall
	for key in functions.keys():
		if not key in functionCall: continue
		remastered = remastered.replace(key, "functions[\""+key+"\"]([X])")

	## now looking for the proper number of combinations that satisfy the condition
	nPassed = 0
	for comb in combinations:
		theObjs = ",".join(["objlist["+str(it)+"]" for it in comb])
		if eval(remastered.replace("X", theObjs)): nPassed += 1

	if nPassed>=req: return True
	return False



def _multi(functions, event, trigObjlists, trigLeg, selectedByOtherLegs, onlyThresholdCut=False, exceptThresholdCut=False):
	## the problem with normal trigger legs is, that if we have more than three such
	## legs acting on the same trigger object (e.g. TkMu), the chronology in how the
	## objects are evaluated and assigned to the individual legs matters; in other
	## words, there may be a combination of objects in the event satisfying the
	## 3 or more leg requirement, but it can be missed by the normal leg since not
	## all of the combinations are evaluated (in order to save time); this function,
	## however, evaluates a multi-leg, thus, builds all combinations in the event
	## and probes if there is at least one combination that meets the requirements

	## collect the trigger objects
	objs  = trigObjlists[trigLeg.obj.name]
	nObjs = len(objs)

	## not enough objects in the event => returning False in form of an empty collection
	if nObjs < trigLeg.num: 
		return [[] for i in range(trigLeg.num)]

	## check if any of the objects is contained in args
	def objsInArgs(theObjs, args):
		for iu in range(len(args)):
			for iv in range(len(args[iu])):
				if any([x in args[iu][iv] for x in theObjs]):
					return True
		return False

	## collecting the objects selected by previously evaluated trigger legs
	args = _getArgs(trigLeg, selectedByOtherLegs)

	## preselection of the objects switching off the part by which the
	## individual sublegs depend on one another; this part is evaluated
	## after the buildling of the combinations
	## this preselection step decreases the number of objects considered
	## for the combinations, thus, can decrease the number of combinations
	## dramatically, which increases computing speed
	goodObjs = []
	for io in range(trigLeg.num):
		goodObjs.append([])
		elms = [goodObjs[i] for i in range(io)]
		for i in range(nObjs):
			if objsInArgs([objs[i]], args): continue
			if not trigLeg.apply(functions, objs[i], io, args, elms, False, onlyThresholdCut=False, exceptThresholdCut=False): continue
			goodObjs[io].append(objs[i])


	## build the combinations of preselected objects
	combinations = []
	nElms = [len(goodObjs[io]) for io in range(trigLeg.num)]
	nVars = int(numpy.prod(nElms))
	for iVar in range(nVars):
		newentry = [(iVar/numpy.prod(nElms[i+1:]))%numpy.prod(nElms[i]) for i in range(len(nElms)-1)] + [iVar%nElms[-1]]
		if any([x!=1 for x in collections.Counter(newentry).values()]): continue
		if newentry in combinations: continue
		combinations.append(newentry)


	## evaluating the rest of the cuts, now considering the full combination
	tried = []
	for comb in combinations:
		theObjs = [goodObjs[io][comb[io]] for io in range(trigLeg.num)]

		#if objsInArgs(theObjs, args): continue
		elms = [[x for x in theObjs[0:io]] for io in range(len(theObjs))]

		## first we check that all members of the combination pass all the
		## cuts in all sublegs, in particular, also the cuts that make the
		## objects depend on one another
		allpass = True
		for io, obj in enumerate(theObjs):
			if not trigLeg.apply(functions, obj, io, args, elms[io], True, onlyThresholdCut=False, exceptThresholdCut=False): 
				allpass = False
				break
		if not allpass: continue

		## the combination now, in principle, is good, but there may be 
		## extra cuts to be evaluated on top of the combination (marked
		## with the $ symbol in the cfg); evaluate them here
		## this is a bit of a problem: suppose you have 3 objects in your
		## combinations; they will be permutations of these three objects
		## while the extra functions do not really care about the 
		## permutations (this is the nature of the extra functions); thus
		## the extra functions could be re-evaluated over and over again
		## only for different permutations; however, you DO want to keep
		## this check here, because the combinations in this loop could
		## contain different set of objects (e.g. there are 5 jets in the
		## event but you build 3-jet combinations); then this is the good
		## place to check the extra functions, but using the "tried" 
		## variable to avoid checking things 
		## it suffices to do all of this once, by the nature of the extra 
		## (functions
		if sum(comb) in tried: continue
		for key,function in trigLeg.extras.iteritems():
			if not function(functions, theObjs):
				allpass = False
		if not allpass: tried.append(sum(comb)); continue


		## if all is good, return the first valid combination
		return [[obj] for obj in theObjs]


	## no appropriate combination has been found, thus we return False 
	## in shape of an empty collection
	return [[] for i in range(trigLeg.num)]






## -------------------------------------------------------------
## -------------------------------------------------------------
## -------------------------------------------------------------
## -------------------------------------------------------------


class TriggerObject:
	## implementation of a trigger object, i.e. the collection of branches
	## read from the tree and used in a trigger leg

	def __init__(self, master, objDef):
		self.master  = master
		self.name    = objDef.name
		self.opts    = objDef.options
		self.load()
	def load(self):
		## parse the cfg and set the variables the trigger object needs 
		basebranch   = self.opts["basebranch"  ] if "basebranch"   in self.opts.keys() else "Muon"
		lengthbranch = self.opts["lengthbranch"] if "lengthbranch" in self.opts.keys() else "n"+basebranch
		separator    = self.opts["separator"   ] if "separator"    in self.opts.keys() else ""
		leadingvar   = self.opts["leadingvar"  ] if "leadingvar"   in self.opts.keys() else "Et"
		leadingop    = self.opts["leadingop"   ] if "leadingop"    in self.opts.keys() else ">"
		isFlat       = self.opts["isFlat"      ] if "isFlat"       in self.opts.keys() else False
		self.onToOff = self.opts["onToOff"     ] if "onToOff"      in self.opts.keys() else 0
		varnames     = []
		varbranch    = []
		for raw in self.opts["variables"]:
			if "=" in raw:
				sr = raw.split("=")
				varnames .append(sr[0])
				varbranch.append(sr[1])
			else:
				varnames .append(raw)
				varbranch.append(raw)
		if not leadingvar in varnames: 
			self.master.vb.warning("Leading variable ("+leadingvar+") does not exist for object "+self.name+"!\nUsing "+varnames[0]+" with operator "+leadingop+" instead. Results may be different than expected!")
			leadingvar = varnames[0]
		self.branches = {}
		for iv, var in enumerate(varnames):
			self.branches[var] = basebranch + separator + varbranch[iv]
		self.lengthbranch = lengthbranch
		self.leadingvar   = leadingvar
		self.leadingop    = leadingop  
		self.isFlat       = isFlat


class TriggerSelection(object):
	## functionality general to any sort of trigger leg (normal or multi)

	def __init__(self):
		self.thresholdCut = ""
		self.thresholdVar = lambda obj: -1
	def buildCut(self, key, cut):
		## the function parses the raw statements of the cfg
		if type(cut)==type(True):
			return key, str(cut), False, []
		if isFloat(cut) or isInt(cut):
			cut = "obj."+self.obj.leadingvar+self.obj.leadingop+str(cut)
			return key, cut, False, []
		if ":" in cut:
			sc = cut.split(":")
			key = sc[0]
			cut = sc[1]
		cut = cut.replace("[","(").replace("]",")")
		uses = list(set([int(i) for i in re.findall(r'\bleg(\d+)\.\b', cut)]))
		elms = list(set([int(i) for i in re.findall(r'\belm(\d+)\.\b', cut)]))
		if not any([x in cut for x in self.obj.opts["variables"]]): 
			self.master.vb.error("Doing something stupid")
		if any([x in cut for x in self.master.functions.keys()]):
			for fname in self.master.functions.keys():
				if fname in cut:
					cut = re.sub(r"\w*(?<![a-zA-Z])(?<!\.)"+fname, "functions[\""+fname+"\"]", cut)
		for var in self.obj.opts["variables"]:
			if var in cut:
				cut = re.sub(r"\w*(?<![a-zA-Z])(?<!\.)"+var, "obj."+var, cut) ## replace only if no letter in front if it (other function?)
		return key, cut, len(elms)>0, uses
	def finalizeCutExpr(self, rawcut):
		values = self.getCutValues(rawcut)
		if len(values)<=1: return rawcut
		var    = self.getCutVar(rawcut)
		ops    = self.getCutOps(rawcut)
		if len(values)==2:
			minval = ops[1][1][0] if ">" in ops[1][0] else ops[0][0][0]
			maxval = ops[0][1][0] if "<" in ops[0][0] else ops[1][0][0]
			betw   = float(minval)<float(maxval)
			eqty   = any(["=" in x for x in ops])
			return ("between" if betw else "outside")+("Eq" if eqty else "")+"("+var+", "+str(minval)+", "+str(maxval)+")"
		return rawcut
	def splitCutString(self, rawcut):
		splitter = None
		if   "<=" in rawcut: splitter = "<="
		elif ">=" in rawcut: splitter = ">="
		elif "<"  in rawcut: splitter = "<"
		elif ">"  in rawcut: splitter = ">"
		elif "!=" in rawcut: splitter = "!="
		elif "==" in rawcut: splitter = "=="
		if splitter:
			return [rawcut.split(splitter), splitter]
		return [[rawcut], None]
	def getCutOps(self, rawcut):
		if type(rawcut)==list:
			return [self.getCutOps(rawcut[i]) for i in range(len(rawcut))]
		last   = [rawcut]
		ops    = []
		save   = []
		result = []
		while True:
			now = []
			for item in last:
				splitted = self.splitCutString(item)
				now.extend (splitted[0])
				if splitted[1]!=None:
					ops.append((splitted[1], filter(lambda x: x[0:len(splitted[1])].isdigit(), splitted[0])))
			if len(now)==len(last): break
			last = now
		return ops
	def getCutValues(self, rawcut):
		## extracting the cut value from the raw cut statement (e.g. 15 out of Et>=15)
		if type(rawcut)==list:
			return [self.getCutValues(rawcut[i]) for i in range(len(rawcut))]
		last = [rawcut]
		while True:
			now = []
			for item in last:
				now.extend(self.splitCutString(item)[0])
			if len(now)==len(last): break
			last = now
		return [float(x) for x in now[1:]] ## we remove the first thing
	def getCutVar(self, rawcut):
		## extracting the cut variable from the raw cut statement (e.g. Et out of Et>=15)
		if type(rawcut)==list:
			return [self.getCutVar(rawcut[i]) for i in range(len(rawcut))]
		last = [rawcut]
		while True:
			now = []
			for item in last:
				now.extend(self.splitCutString(item)[0])
			if len(now)==len(last): break
			last = now
		return now[0] ## we keep only the first
	def setCutValue(self, rawcut, newthreshold):
		## replacing the cut value in the raw cut statement (e.g. Et>=15 to Et>=20)
		if type(rawcut)==list and len(rawcut)==len(newthreshold):
			return [self.setCutValue(rawcut[i], newthreshold[i]) for i in range(len(rawcut))]
		if not any([x in rawcut for x in ["<",">","="]]): return float(newthreshold)
		if "<=" in rawcut: return rawcut.split("<=")[0]+"<="+str(newthreshold)
		if ">=" in rawcut: return rawcut.split(">=")[0]+">="+str(newthreshold)
		if "<"  in rawcut: return rawcut.split("<" )[0]+"<" +str(newthreshold)
		if ">"  in rawcut: return rawcut.split(">" )[0]+">" +str(newthreshold)
		if "!=" in rawcut: return rawcut.split("!=")[0]+"!="+str(newthreshold)
		if "==" in rawcut: return rawcut.split("==")[0]+"=="+str(newthreshold)
		return float(rawcut)
	def setThresholdCut(self, name):
		## determining what is the threshold variable for this leg (probably it's the Et)
		## and compiling a function that will return the objects value for this variable
		## when applied to an arbitrary object
		self.thresholdCut = name
		cutvar            = self.getCutVar(self.cuts[name])
		self.thresholdVar = eval("lambda obj: "+cutvar if cutvar else self.obj.leadingvar)



class TriggerLeg(TriggerSelection):
	## implementation of a normal trigger leg, i.e. a single selection applied to an object

	def __init__(self, master, obj, cuts, idx):
		super(TriggerLeg,self).__init__()
		self.master   = master
		self.obj      = obj
		self.idx      = int(idx)
		self.raw      = " ".join([str(c) for c in cuts])
		self.uses     = []
		self.buildDef(cuts)
		self.characterize()
	def getThreshold(self, obj):
		## return the value of an object for the variable that is the threshold variable 
		## for this trigger leg (it's probably the Et, but can be changed)
		return self.thresholdVar(obj)
	def apply(self, functions, obj, args, onlyThresholdCut=False, exceptThresholdCut=False):
		## evaluating the leg on an object (obj)

		## loop over all individual cuts in the leg
		for name,function in self.functs.iteritems():
			isThresholdCut = (name==self.thresholdCut)
			if (onlyThresholdCut and not isThresholdCut) or (exceptThresholdCut and isThresholdCut): continue
			iu = int(name.replace("cut",""))

			## the leg depends on a previous leg but the corresponding objects are not given
			if len(self.uses[iu])>0 and len(args[iu])==0: 
				return False

			## without dependency on previous legs
			if len(args[iu])==0:
				if not function(functions, obj): return False

			## with dependency on previous legs
			else:
				## args[iu] contains all variations of the list of arguments
				## only one variation needs to work to trigger the event
				## N.B. the object selected here, if in a different collection
				## than previous legs (e.g. leg1=TkEG, leg2=EG), can still
				## overlap with the object selected at any of the previous legs
				## but this is no problem since at least one combination of objects
				## needs to work; when probing further legs (e.g. leg3 depending
				## on leg1 and leg2), the proper combination to probe is among
				## the ones that is probed
				anytrue = False
				for arg in args[iu]:
					if len(self.uses[iu])!=len(arg): continue
					if function(functions, obj, *arg): anytrue=True
				if not anytrue: return False
		return True
	def buildDef(self, cuts):
		## building the trigger from the cfg definition
		self.cuts = {}
		for ic,cut in enumerate(cuts):
			key = "cut"+str(ic)
			self.uses.append([])
			key, cut, elms, self.uses[ic] = self.buildCut(key, cut)
			self.cuts[key] = cut
		## compiling the individual cuts (in self.functs) that are applied/evaluated per object
		self.functs  = {}
		for key,cut in self.cuts.iteritems():
			iu = int(key.replace("cut",""))
			toUse = self.uses[iu]
			self.functs[key] = eval("lambda functions, obj"+(", "+", ".join("leg"+str(i) for i in toUse) if len(toUse)>0 else "")+": "+cut)
	def characterize(self):
		## find and store the id of the virtual leg
		self.legId = [self.master.getLegId(["leg", self.obj.name] + [c for k,c in self.cuts.iteritems()])]


class TriggerMulti(TriggerSelection):
	## implementation of a multi-leg, i.e. a leg with num>1 "sublegs"
	## each depending on the same object obj; recommended to use if
	## num>2 (if num=2 better use two individual normal legs) since
	## it explores all possible combinations of the objects

	def __init__(self, master, obj, num, cuts, idx):
		super(TriggerMulti,self).__init__()
		self.master   = master
		self.obj      = obj
		self.num      = int(num)
		self.idx      = int(idx)
		self.raw      = " ".join([str(c) for c in cuts])
		self.uses     = []
		self.elms     = []
		self.buildDef(cuts)
		self.characterize()
	def getThreshold(self, obj):
		## return the value of an object for the variable that is the threshold variable 
		## for this trigger leg (it's probably the Et, but can be changed)
		return self.thresholdVar(obj)
	def apply(self, functions, obj, io, args, elms, onlyWithElms=True, onlyThresholdCut=False, exceptThresholdCut=False):
		## evaluating the io'th part of the multi-leg on an object (obj)

		## loop over all individual cuts in the leg
		for name,function in self.functs.iteritems():
			iu = int(name.replace("cut",""))
			isThresholdCut = (name==self.thresholdCut)
			if (onlyThresholdCut and not isThresholdCut) or (exceptThresholdCut and isThresholdCut): continue
			if (onlyWithElms and not self.elms[iu]) or (not onlyWithElms and self.elms[iu]): continue

			## the leg depends on a previous leg but the corresponding objects are not given
			if len(self.uses[iu])>0 and len(args[iu])==0: 
				return False 

			## without dependency on previous legs
			if len(args[iu])==0:
				if len(elms)==0:
					if not function[io](functions, obj): return False
				else:
					if not function[io](functions, obj, *elms): return False

			## with dependency on previous legs
			else:
				## args[iu] contains all variations of the list of arguments
				## only one variation needs to work to trigger the event
				## N.B. the object selected here, if in a different collection
				## than previous legs (e.g. leg1=TkEG, leg2=EG), can still
				## overlap with the object selected at any of the previous legs
				## but this is no problem since at least one combination of objects
				## needs to work; when probing further legs (e.g. leg3 depending
				## on leg1 and leg2), the proper combination to probe is among
				## the ones that is probed
				anytrue = False
				for arg in args[iu]:
					if len(self.uses[iu])!=len(arg): continue
					myargs = arg + elms
					if function[io](functions, obj, *myargs): anytrue=True
				if not anytrue: return False
		return True
	def buildDef(self, cuts):
		## building the trigger from the cfg definition
		self.cuts = {}
		self.exts = {}
		for ic,cut in enumerate(cuts):
			key = "cut"+str(ic)
			self.uses.append([])
			self.elms.append([])
			if type(cut)==list:
				self.cuts[key] = [True for i in range(self.num)]
				for ii, c in enumerate(cut):
					key, cutstring, self.elms[ic], self.uses[ic] = self.buildCut(key, c)
					self.cuts[key][ii] = self.finalizeCutExpr(cutstring)
			else:
				isExtra = False
				if type(cut)==type("bla") and cut[0:1]=="$":
					isExtra = True
					ext = cut[cut.find("[")+1:cut.find("]")].split(":")
					req = int(ext[0]) ## how many combinations
					num = int(ext[1]) ## how many objects for each combination
					cut = self.finalizeCutExpr(cut[cut.find("]")+1:])
					self.exts[key] = (cut, num, req)
					continue	
				self.cuts[key] = [True for i in range(self.num)]
				key, cutstring, self.elms[ic], self.uses[ic] = self.buildCut(key, cut)
				for ii in range(self.num):
					self.cuts[key][ii] = self.finalizeCutExpr(cutstring)
		## compiling the individual cuts (in self.functs) and the extra cuts (in self.extras)
		## the former are used per object in order to determine the proper combination of
		## objects; the latter are used on the entire combination
		self.functs = {}
		self.extras = {}
		for key,cut in self.cuts.iteritems():
			self.functs[key] = [{} for i in range(self.num)]
			iu = int(key.replace("cut",""))
			toUseL = self.uses[iu]
			for ic in range(self.num):
				astring = ", ".join("leg"+str(i) for i in toUseL       ) if len(toUseL)>0 else ""
				estring = ", ".join("elm"+str(i) for i in range(1,ic+1)) if ic>0          else ""
				full = ", "+astring+", "+estring if astring!="" and estring!="" else ", "+astring if astring!="" else ", "+estring if estring!="" else ""
				self.functs[key][ic] = eval("lambda functions, obj"+full+": "+cut[ic])
		objs = ",".join(["obj"+str(i) for i in range(self.num)])
		for key,cut in self.exts.iteritems():
			self.extras[key] = eval("lambda functions, objList: _multiPermute(functions, \""+cut[0]+"\", "+str(cut[1])+", "+str(cut[2])+", objList)")
	def characterize(self):
		## find and store the id of the virtual multi-leg
		self.legId = []
		for io in range(self.num):
			self.legId.append(self.master.getLegId(["multi"+str(io), self.obj.name] + [c[io] for k,c in self.cuts.iteritems()]))



class Trigger:
	## implementation of the entire trigger with all its legs
	## the individual entities and functions are built from the
	## cfg definition which is encoded in the trigDef variable

	def __init__(self, master, trigDef):
		self.master   = master
		self.name     = trigDef.name
		self.opts     = trigDef.options
		self.parent   = None
		self.trigVar  = {}
		self.load()
		self.characterize()
	def apply(self, event, trigObjlists, onlyThresholds=False, exceptThresholds=False):
		## applying the trigger on an event with the corresponding trigger objects
		## given in the trigObjlists variable
		legs = []
		for leg in self.legs.values():
			legs.append(_leg(self.master.functions, event, trigObjlists, leg, legs))
		for mul in self.multis.values():
			multi = _multi(self.master.functions, event, trigObjlists, mul, legs)
			legs.extend(multi)
		## the legs variable contains a list of objects passing the individual legs
		## if any of these lists is empty, it means, the event failed the corresponding leg
		if any([len(l)==0 for l in legs]): return False
		return True 
	def characterize(self):
		## find and store the virtual id of the trigger
		## this relies on the virtual id of the individual legs
		legIds = []
		for leg in self.legs.values():
			legIds.extend(leg.legId)
		for leg in self.multis.values():
			legIds.extend(leg.legId)
## CH: check this look up of the trigger id again, especially for the variation case!!
		self.triggerId = self.master.getTrigId(["single"] + legIds)
	def getCutValues(self, legname, cutname):
		if not legname in self.legnames: return None
		collection = self.multis if "multi" in legname else self.legs
		if not cutname in collection[legname].cuts.keys(): return None
		rawExpr = collection[legname].cuts[cutname]
		return collection[legname].getCutValues(rawExpr)
	def getLegNameByIdx(self, legidx):
		if legidx>=len(self.legnames): return None
		return self.legnames[legidx]
	def getLegByName(self, legname):
		if not legname in self.legnames: return None
		return self.legs[legname] if "leg" in legname else self.multis[legname]
	def makeVar(self, varname, varleg, varcut, varvalue):
		## create a copy of the trigger with a varied cut
		## i.e. the cuts for ALL legs for this trigger are varied
		## note also that all sublegs need to be varied accordingly for every
		## multi-leg that the trigger might contain

		newDef    = CfgObject(self.name+"_"+varname,"")
		varvalue  = float(varvalue)
		legToVary = self.legnames[0] if varleg=="any" else varleg
		theLeg    = self.multis[legToVary] if "multi" in legToVary else self.legs[legToVary]
		add       = 2 if "multi" in legToVary else 1 ## entry 0 = object, entry 1 = first cut, +1 in multi
		cutIdx    = int(varcut.replace("cut",""))-1
		cutToVary = "cut"+str(cutIdx)
		#cutToVary = int(varcut.replace("cut",""))+add

		## when multiple legs are present, varleg and varvalue will only affect one of the legs
		## thus need to keep the ratio of these cuts fixed and vary all other legs accordingly !!!
		oldvalues = [self.getCutValues(ln, cutToVary) for ln in self.legnames]
		newvalues = oldvalues[:]
		refs      = oldvalues[self.legnames.index(legToVary)]
		varvalues = [varvalue]
		if type(refs)!=list: refs = [refs]
		if "multi" in legToVary and len(refs)>1: 
			## if only one value is given for a multi-leg, it means that it is the leading value, 
			## and all other legs must be scaled up too accordingly 
			for ie in range(1,len(refs)):
				varvalues.append(refs[ie]*(varvalue/refs[0]))
		scales = [varvalues[i]/refs[i] for i in range(len(varvalues))]
		if sum(scales)==len(scales): 
			self.varValue = varvalue
			return self
		for iv,oldvalue in enumerate(oldvalues):
			if type(oldvalue)==list:
				newvalues[iv] = []
				for ie,elm in enumerate(oldvalue):
					newvalues[iv].append(theLeg.setCutValue(cutToVary, scales[ie]*elm))
				continue
			newvalues[iv] = theLeg.setCutValue(cutToVary, scales[0]*oldvalue) 
			## this is necessary in case the cut you want to vary contains an expression, not just a float or int value
		for key in self.opts.keys():
			if not key in self.legnames: 
				newDef.options[key] = self.opts[key]
				continue
			idx = self.legnames.index(key)
			newleg              = self.opts[key]
			newleg[cutIdx+add]  = newvalues[idx]
			newDef.options[key] = newleg
		newTrig = Trigger(self.master, newDef)
		newTrig.parent   = self
		newTrig.varName  = varname
		newTrig.varLeg   = varleg  
		newTrig.varCut   = varcut  
		newTrig.varValue = varvalue
		self.trigVar[varname] = newTrig
		return newTrig
	def load(self):
		## prepare the functions and all other constraints
		## builds list of trigger legs and event selection, each of them functions to be 
		## executed the apply(evt, trigObjlists) then runs all of them with trigObjlists the 
		## actual Objectlists instantiated for a given event
		## for a leg, it takes a _leg(evt, trigObjlists, legDef)
		## for a evt, it takes a _evt(evt, trigObjlists, legs  , evtDef)
		objnames         = [o.name for o in self.master.menu.objects]
		self.legs        = {}
		self.multis      = {}
		self.ranges      = {}
		self.legnames    = []
		self.legentities = []
		for key in self.opts.keys():
			## a normal leg
			if "leg" in key:
				args = self.opts[key]
				if not args[0] in objnames:
					self.master.vb.warning("Cannot find trigger object "+args[0]+" for "+key+" of path "+self.name+"\nSkipping this leg! Results may be different than expected!")
					continue
				obj = self.master.menu.objects[objnames.index(args[0])]
				self.legs[key] = TriggerLeg(self.master, obj, args[1:], int(key[3:]))
				self.legnames   .append(key)
				self.legentities.append(self.legs[key])

			## a multi-leg
			if "multi" in key:
				args = self.opts[key]
				if not args[0] in objnames:
					self.master.vb.warning("Cannot find trigger object "+args[0]+" for "+key+" of path "+self.name+"\nSkipping this leg! Results may be different than expected!")
					continue
				obj = self.master.menu.objects[objnames.index(args[0])]
				self.multis[key] = TriggerMulti(self.master, obj, args[1], args[2:], int(key[5:]))
				self.legnames   .append(key)
				self.legentities.append(self.multis[key])

			if "evt" in key:
				pass ## not implemented currently

			## a range to vary the threshold of a previously defined trigger leg
			if "range" in key:
				args = self.opts[key]
				if len(args)<3:
					self.master.vb.warning("Not enough arguments for "+key+" of path "+self.name+"\nSkipping this range!")
					continue
				if not args[0] in self.legs.keys() and not args[0] in self.multis.keys():
					self.master.vb.warning("Cannot find leg "+args[0]+" for "+key+" of path "+self.name+"\nSkipping this range!")
					continue
				self.ranges[key] = args
	def setThresholdCuts(self):
		varCut = "cut"+str(self.master.cfg.variable["varCut"])
		for leg   in self.legs  .values():
			leg  .setThresholdCut(varCut)
		for multi in self.multis.values():
			multi.setThresholdCut(varCut)
	def setThresholdValues(self, legsresults):
		self.thresholdValues = []
		for i,legresult in enumerate(legsresults):
			self.thresholdValues.append([])
			theLeg = self.legentities[i]
			if len(legresult)==0: self.thresholdValues[i] = [-1]; continue
			for legobj in legresult:
				self.thresholdValues.append(theLeg.getThreshold(legobj))
	def getThresholdValues(self, trig):
		return self.thresholdValues


	
