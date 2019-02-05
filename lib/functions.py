import ROOT, re, os, datetime, math
import inspect, pandas

def cmd(base):
	#print base
	os.system(base)

def mkdir(path):
	if os.path.isdir(path): return
	cmd("mkdir -p "+path)

def cleandir(path):
	if not os.path.isdir(path): return
	cmd("rm -rf "+path+"/*")

def rmFile(path):
	if not os.path.exists(path): return
	cmd("rm -rf "+path)

def isFloat(string):
	if type(string)==type(1.0)  : return True
	if type(string)!=type("bla"): return False
	if not "." in string        : return False
	try:
		return float(string)
	except (ValueError, TypeError) as e:  
		return False

def isInt(string):
	if type(string)==type(1)    : return True
	if type(string)!=type("bla"): return False
	if "." in string            : return False
	try:
		return float(string)
	except (ValueError, TypeError) as e:  
		return False

def getRootObj(tfile, objDef):
	if not "/" in objDef: return tfile.Get(objDef)
	so = objDef.split("/")
	prev = tfile.Get(so[0])
	for sk in so[1:]:
		if prev.ClassName()=="TDirectoryFile" and sk==so[-1]:
			## FIXME: could be improved
			ttree = ROOT.TTree()
			prev.GetObject(sk, ttree)
			return ttree
		prev = prev.Get(sk) if "Get" in prev.__dict__ and prev.Get(sk) else prev.GetPrimitive(sk)
	return prev

def getType(variable, isList=False, isDict=False):
	r = re.compile(r'(?:[^,(]|\([^)]*\))+')
	if variable[0:1]=="(" and variable[-1:]==")":
		variable = variable[1:-1]
		isList   = True
	if isDict:
		cache = [ss.strip().split(":") for ss in r.findall(variable)]
		return (dict, {ss[0].strip():getType(ss[1].strip())[1] for ss in cache})
	if isList:
		return (list, [getType(ss.strip())[1] for ss in r.findall(variable)])
		#return (list, [getType(ss.strip())[1] for ss in variable.split(",")])
	if variable=="True": 
		return (bool, True)
	if variable=="False": 
		return (bool, False)
	if type(variable)==type(True):
		return (bool, bool(variable))
	if isInt(variable):
		return (int, int(variable))
	if isFloat(variable):
		return (float, float(variable))
	return (basestring, variable)

def setType(variable, isList=False, isDict=False):
	tple = getType(variable, isList, isDict)
	return tple[1]

def shortstamp(readable = True):
	if readable:
		return datetime.datetime.now().strftime("%Y-%m-%d")
	return datetime.datetime.now().strftime("%y%m%d")

def timestamp(readable = True):
	if readable:
		return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	return datetime.datetime.now().strftime("%y%m%d%H%M%S%f")

def absPath(relative, absolute):
	if relative[0:1]=="/": return relative.rstrip("/")
	return absolute.rstrip("/")+"/"+relative.rstrip("/")

def mergeLines(lines):
	newlines = []
	buff     = ""
	for line in lines:
		if line.strip()[0:1]=="#": continue
		if line[-1:]=="\\": buff+=line[0:-1]; continue
		newlines.append(buff+line)
		if line[-1:]!="\\": buff=""
	return newlines

def safeDiv(numerator, denominator):
	if denominator == 0: return 0
	return float(numerator)/denominator

def divLists(numerators, denominators):
	if len(numerators)!=range(denominators): return []
	quotients = []
	for i in range(len(numerators)):
		quotients.append(safeDiv(numerators[i], denominators[i]))
	return quotients

def readObjsFromRootFile(rfile):
	contents = []
	rfile.cd()
	for key in rfile.GetListOfKeys():
		obj = rfile.Get(key.GetName())
		if obj and "SetDirectory" in [m[0] for m in inspect.getmembers(obj, predicate=inspect.ismethod)]:
			obj.SetDirectory(0)
		contents.append(obj)
	return contents

def readCache(path, idx="1", selection = lambda x: True):
	collection = {}
	if not os.path.exists(path): return collection
	f = open(path, "r")
	for line in [l.rstrip("\n").strip() for l in f.readlines()]:
		if len(line) == 0: continue
		if line[0:1] == "#": continue
		fields = [sl.strip() for sl in line.split("::")]
		if len(fields)<2        : continue
		if not selection(fields): continue
		collection[fields[0]] = fields[int(idx.isdigit())] if idx.isdigit() else eval("fields["+idx+"]")
	f.close()
	return collection

def writeCache(path, collection, asList = True):
	if len(collection.keys())==0: return
	rmFile(path)
	f = open(path,"w")
	for k,v in collection.iteritems():
		if asList: f.write(k+" :: "+" :: ".join(v)+"\n")
		else     : f.write(k+" :: "+v+"\n")
	f.close()

def updateTriggers(fullSet, varTrigger):
	theSet = []
	for trig in fullSet:
		if trig == varTrigger.parent: theSet.append(varTrigger); continue
		theSet.append(trig)
	return theSet

def findIntersectionIdx(collection, element):
	## find the idx of the element in the collection that has value just below (above) 
	## that of the seeked element, i.e. idx+1 will have a value above (below) it
	previous = 0
	if element>collection[0                ]: return -1
	if element<collection[len(collection)-1]: return len(collection)
	for i,x in enumerate(collection):
		this = x-element ##>0 if rate is smaller than x, <0 if rate is larger than x
		if previous*this<0: ## change sign
			return i-1
			break
		previous = this
	return -1

def getPoissonErrors(values, errors=[]):
	q  = (1-0.6827)/2. # quantile
	er = [math.sqrt(x) for x in values] if len(errors)==0 else errors
	up = []
	dn = []
	for i,value in enumerate(values):
		N  = value
		dN = er[i]
		if N==0: up.append(0); dn.append(0); continue
		scale = 1
		if N>0 and dN>0 and abs(dN**2/N-1)>1e-4:
			scale = (dN**2/N)
			N     = (N/dN)**2
		up.append(scale*((N-ROOT.ROOT.Math.chisquared_quantile_c(1-q,2*N)/2.) if N>0 else 0))
		dn.append(scale*(ROOT.ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2.-N))
	return up, dn

def binomial(n, k):
	if n == k: return 1
	if k == 1: return n
	if k >  n: return 0
	return math.factorial(n) // (math.factorial(n-k) * math.factorial(k))

def between(x, xmin, xmax):
	return xmin<x and x<xmax

def betweenEq(x, xmin, xmax):
	return xmin<=x and x<=xmax

def outside(x, xmin, xmax):
	return x<xmax or x>xmin

def outsideEq(x, xmin, xmax):
	return x<=xmax or x>=xmin

def concatDFs(dfs, sels, fields, how):
	if len(dfs)==1: return dfs[0], sels[0]
	## first the big one without selection 
	df = reduce(lambda left, right: pandas.merge(left, right, on=fields, how=how), dfs)
	df.dropna(inplace=True)
	## then the small one with all the selections
	df2 = reduce(lambda left, right: pandas.merge(left, right, on=fields, how=how), sels)
	df2.dropna(inplace=True)
	return df, df2

def subsetDF(df, sel):
	dfsel = df.loc[df['ievt'].isin(sel)]
	#dfsel.dropna(inplace=True)
	return dfsel

def intersectDFs(dfs, sels, fields):
	return concatDFs(dfs, sels, fields, "inner")

def mergeDFs(dfs, sels, fields):
	if len(dfs)==1: return dfs[0], sels[0]
	## first the big one without selection 
	df = pandas.concat(dfs, sort=False)
	df.drop_duplicates(subset="ievt")
	df.reset_index(drop=True)
	## then the small one with all the selections
	df2 = pandas.concat(sels, sort=False)
	df2.drop_duplicates(subset="ievt")
	df2.reset_index(drop=True)
	return df, df2

#	return concatDFs(dfs, sels, fields, "inner")


