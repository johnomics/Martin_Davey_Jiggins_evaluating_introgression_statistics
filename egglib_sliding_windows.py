#!/usr/bin/env python

# egglib_sliding_windows.py
# calculkate ABBA BABA stats, dxy and S for sliding windows in genomic data

# Written for "Evaluating statistics for the identification of introgressed loci"
# by Simon H. Martin, John W. Davey and Chris D. Jiggins
# John Davey:   jd626@cam.ac.uk
# Simon Martin: shm45@cam.ac.uk
# November-December 2013



import sys

import egglib

def getOptionValue(option): # needs sys
  if option in sys.argv:
    optionPos = sys.argv.index(option)
    optionValue = sys.argv[optionPos + 1]
    return optionValue
  else:
    print >> sys.stderr, "\nWarning, option", option, "not_specified.\n"


def get_intv(string,borders = "()",inc = False):
  if len(borders) != 2:
    print "WARNING: borders must contain two characters"
  starts = []
  ends = []
  output = []
  for x in range(len(string)):
    if string[x] == borders[0]:
      starts.append(x)
    if string[x] == borders[1]:
      ends.append(x+1)
  if len(starts) <= len(ends):
    for n in range(len(starts)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  else:
    for n in range(len(ends)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  return output

def median(numbers):
  numbers.sort()
  if len(numbers) % 2 == 1:
    return numbers[(len(numbers)+1)/2-1]
  else:
    lower = numbers[len(numbers)/2-1]
    upper = numbers[len(numbers)/2]
    return (float(lower + upper)) / 2

def haplo(calls):
  output = []
  for call in calls:
    if call in "ACGTN":
      output.append(call)
      output.append(call)
    elif call == "K":
      output.append("G")
      output.append("T")
    elif call == "M":
      output.append("A")
      output.append("C")
    elif call == "R":
      output.append("A")
      output.append("G")
    elif call == "S":
      output.append("C")
      output.append("G")
    elif call == "W":
      output.append("A")
      output.append("T")
    elif call == "Y":
      output.append("C")
      output.append("T")
    else:
      print "WARNING", call, "is not recognised as a valid base or ambiguous base"
      output.append("N")
      output.append("N")
  return output

def mean(numbers):
  numbers = [float(n) for n in numbers if n != "NA" and n != None]
  numSum = sum(numbers)
  if len(numbers) >= 1:
    return float(numSum)/len(numbers)
  else:
    return "NA"

def mid(numbers):
  numbers = [float(n) for n in numbers if n != "NA" and n != None]
  if len(numbers) >= 1:
    return (numbers[0] + numbers[-1])/2
  else:
    return None

def AlignByGroupNumber(align,groupNumber):
  newAlign = align.slice(0,0)
  for seqNumber in range(len(align)):
    if align[seqNumber][2] == groupNumber:
      newAlign.addSequences([align[seqNumber]])
  return newAlign

def AlignByGroupNumbers(align,groupNumbers):
  newAlign = align.slice(0,0)
  for seqNumber in range(len(align)):
    if align[seqNumber][2] in groupNumbers:
      newAlign.addSequences([align[seqNumber]])
  return newAlign

def mostCommon(things):
  output = []
  counts = []
  uniqueThings = unique(things)
  for thing in uniqueThings:
    counts.append(things.count(thing))
  maxCount = max(counts)
  for n in range(len(counts)):
    if counts[n] == maxCount:
      output.append(uniqueThings[n])
  return output

def exclude(things, x):
  return [i for i in things if i != x]

def unique(things):
  output = list(set(things))
  output.sort()
  return output


def dxy(align): # "align" if the egglib alignment object, this consistes of sequences, sequence names and "groups". If the object contains two groups, the function will consider only the first two.
    
    # retrieve group names from the alignment object
    pops = align.groups().keys()
    
    # retrieve all the positions of sequences in group 1
    P1 = [i for i in range(len(align)) if align.group(i)==pops[0]]
    
    # retrieve all the positions of sequences in group 2
    P2 = [i for i in range(len(align)) if align.group(i)==pops[1]]
    
    pairwiseSum = 0 #total of pairwise Pis
    totalPairs = 0 #haplotype pairs considered
    
    for i in P1: # for each sequence in pop1...
      for j in P2: #for sequence in pop2...
        seqA = align[i][1]
        seqB = align[j][1]
        zippedSeqs = zip(seqA,seqB)
        diffs = sum(sA != sB for sA, sB in zippedSeqs if sA != "N" and sB != "N")
        #sites = sum(sA != "N" and sB != "N" for sA, sB in zippedSeqs)
        sites = len([site for site in zippedSeqs if site[0] != "N" and site[1] != "N"])
                  
    #after considering all positions for each pair of haplotypes, return the average pairwise pi
    return 1.0 * diffs/sites


def px(align):
    
    pairwiseSum = 0 #total of pairwise Pis
    totalPairs = 0 #haplotype pairs considered
    
    for i in range(len(align) - 1): # for each sequence except the last one...
      for j in range(i + 1,len(align)): #for each of the remaining sequences from sequence i + 1 to the end of the alignment... 
        seqA = align[i][1]
        seqB = align[j][1]
        zippedSeqs = zip(seqA,seqB)
        diffs = sum(sA != sB for sA, sB in zippedSeqs if sA != "N" and sB != "N")
        #sites = sum(sA != "N" and sB != "N" for sA, sB in zippedSeqs)
        sites = len([site for site in zippedSeqs if site[0] != "N" and site[1] != "N"])
        #now add this pairwise pi to the total and add 1 to the number of pairs considered
        pairwiseSum += 1.0*diffs/sites
        totalPairs += 1
                  
    #after considering all positions for each pair of haplotypes, return the average pairwise pi
    return pairwiseSum/totalPairs


def colFreqs(align, columnNumber):
  bases = align.column(columnNumber)
  Acount = float(bases.count("A"))
  Ccount = float(bases.count("C"))
  Gcount = float(bases.count("G"))
  Tcount = float(bases.count("T"))
  total = Acount + Ccount + Gcount + Tcount
  if total > 0:
    output = {}
    output["A"] = Acount/total
    output["C"] = Ccount/total
    output["G"] = Gcount/total
    output["T"] = Tcount/total
  else:
    output = {"A":"NA", "C":"NA", "G":"NA", "T":"NA"}
  return output

def colBaseCounts(align, columnNumber):
  output = {}
  bases = align.column(columnNumber)
  Acount = float(bases.count("A"))
  Ccount = float(bases.count("C"))
  Gcount = float(bases.count("G"))
  Tcount = float(bases.count("T"))
  output["A"] = Acount
  output["C"] = Ccount
  output["G"] = Gcount
  output["T"] = Tcount
  return output



#version using frequencies as in Durand et al 2011 eqn 2.
def ABBABABA(align, P1, P2, P3, P4):
  p1Align = AlignByGroupNumber(align,P1)
  p2Align = AlignByGroupNumber(align,P2)
  p3Align = AlignByGroupNumber(align,P3)
  p4Align = AlignByGroupNumber(align,P4)
  ABBAsum = 0.0
  BABAsum = 0.0
  maxABBAsum = 0.0
  maxBABAsum = 0.0
  maxABBAsumB = 0.0
  maxBABAsumB = 0.0
  #get derived frequencies for all biallelic siites
  for i in align.polymorphism(minimumExploitableData = 0)["siteIndices"]:
    #skip this site if not biallelic
    bases = exclude(align.column(i), "N")
    alleles = unique(bases)
    if len(alleles) != 2: continue
    #get derived state
    #if the outgroup is fixed, then that is the ancestral state - otherwise the anc state is the most common allele overall
    p4Alleles = unique(exclude(p4Align.column(i), "N"))
    if len(p4Alleles) == 1:
      derived = [a for a in alleles if a != p4Alleles[0]][0]
    else:
      derived = [a for a in alleles if a != mostCommon(bases)[0]][0]
    # get frequencies for wach pop
    p1Freq = colFreqs(p1Align, i)[derived]
    p2Freq = colFreqs(p2Align, i)[derived]
    p3Freq = colFreqs(p3Align, i)[derived]
    p4Freq = colFreqs(p4Align, i)[derived]
    # get weigtings for ABBAs and BABAs
    try: # this was added to ignore crashes when there is missing data for a population at a site - we just ignore these sites
      ABBAsum += (1 - p1Freq) * p2Freq * p3Freq * (1 - p4Freq)
      BABAsum += p1Freq * (1 - p2Freq) * p3Freq * (1 - p4Freq)
      maxABBAsum += (1 - p1Freq) * p3Freq * p3Freq * (1 - p4Freq)
      maxBABAsum += p1Freq * (1 - p3Freq) * p3Freq * (1 - p4Freq)
      if p3Freq >= p2Freq:
        maxABBAsumB += (1 - p1Freq) * p3Freq * p3Freq * (1 - p4Freq)
        maxBABAsumB += p1Freq * (1 - p3Freq) * p3Freq * (1 - p4Freq)
      else:
        maxABBAsumB += (1 - p1Freq) * p2Freq * p2Freq * (1 - p4Freq)
        maxBABAsumB += p1Freq * (1 - p2Freq) * p2Freq * (1 - p4Freq)
    except:
      continue
  #calculate D, f and fb
  output = {}
  try:
    output["D"] = (ABBAsum - BABAsum) / (ABBAsum + BABAsum)
  except:
    output["D"] = "NA"
  try:
    output["f"] = (ABBAsum - BABAsum) / (maxABBAsum - maxBABAsum)
  except:
    output["f"] = "NA"
  try:
    output["mf"] = (ABBAsum - BABAsum) / (maxABBAsumB - maxBABAsumB)
  except:
    output["mf"] = "NA"
  output["ABBA"] = ABBAsum
  output["BABA"] = BABAsum
  output["maxABBA"] = maxABBAsum
  output["maxBABA"] = maxBABAsum
  output["maxABBA_B"] = maxABBAsumB
  output["maxBABA_B"] = maxBABAsumB
  
  return output

#***************************************************************************************************************


if "--stop-at" in sys.argv:
  stopAt = True
  stopVal = int(getOptionValue("--stop-at"))
else:
  stopAt = False

if "--test" in sys.argv:
  test = True
else:
  test = False

if "--verbose" in sys.argv:
  verbose = True
else:
  verbose = False

if "--report" in sys.argv:
  report = int(getOptionValue("--report"))
else:
  report = 100

if "--haplotypes" in sys.argv:
  haplotypes = True
  genotypes = False
  print "\nCalls are: haplotypes with phase\n"
else:
  haplotypes = False
  genotypes = True
  print "\nCalls are: genotypes without phase\n"


if "-i" in sys.argv:
  fileName = getOptionValue("-i")
else:
  print "\nplease specify input file name using -i <file_name> \n"
  sys.exit()

file = open(fileName, "rU")
#define names from header line (file must have a header)
line = file.readline()
names = line.split()
line= file.readline()


if "-p" in sys.argv:
  popStrings = getOptionValue("-p")
else:
  print "\nplease specify populations using -p\n"
  sys.exit()



popNames = []
indNames = []
#for each pattern, store the name, the set of lists, the maximum Ns and the maximum mismatches
for popString in popStrings.strip("\"").split(";"):
  currentPop = popString.split("[")[0]
  popNames.append(currentPop)
  vars()[currentPop] = get_intv(popString,"[]")[0].split(",")
  for indName in vars()[currentPop]:
    if indName in names:
      if indName not in indNames:
        indNames.append(indName)
    else:
      print "individual " + indName + "not found in header line."
      sys.exit()

if "-O" in sys.argv:
  includeOutGroup = True
  outGroupString = getOptionValue("-O").strip("\"")
  outGroup = outGroupString.split("[")[0]
  vars()[outGroup] = get_intv(outGroupString,"[]")[0].split(",")
  for indName in vars()[outGroup]:
    if indName in names:
      indNames.append(indName)
    else:
      print "individual " + indName + "not found in header line."
      sys.exit()
else:
  includeOutGroup = False
  

if test or verbose:
  print "\nPopulations:\n"
  for popName in popNames:
    print popName
    print vars()[popName]
    print "\n"
  if includeOutGroup:
    print "\nOut-Group:\n"
    print outGroup
    print vars()[outGroup]
    print "\n"

#set up a variable that reports the ploidy for each individual

ploidy = {}

if "--ploidy" in sys.argv:
  ploidyNumbers = getOptionValue("--ploidy").strip("\"").split(",")
  ploidyNumbers = [int(n) for n in ploidyNumbers if n == "1" or n == "2"]
  if len(ploidyNumbers) == len(indNames):
    print "\nPloidy is as follows:\n"
    for x in range(len(indNames)):
      ploidy[indNames[x]] = ploidyNumbers[x]
      print indNames[x], ploidyNumbers[x]
  else:
    print "\nSpecify ploidy for each individual as 1 or 2, separated by commas\n"
    sys.exit()
else:
  #if ploidy is not specified, assume diploid
  for indName in indNames:
      ploidy[indName] = 2


if "-o" in sys.argv:
  outName = getOptionValue("-o")
else:
  print "\nplease specify output file name using -o <file_name> \n"
  sys.exit()


if "-w" in sys.argv:
  windSize = int(getOptionValue("-w"))
else:
  print "\nplease specify window size using -w \n"
  sys.exit()

if "-s" in sys.argv:
  slide = int(getOptionValue("-s"))
else:
  print "\nplease specify slide length using -s \n"
  sys.exit()

if "-m" in sys.argv:
  minSites = int(getOptionValue("-m"))
else:
  print "\nplease specify the minimum number of sites per window using -m\n"
  sys.exit()

if "-S" in sys.argv:
  scafsOfInterest = getOptionValue("-S").strip("\"").split(",")
  if test or verbose:
      print "scaffolds to analyse:", scafsOfInterest
  allScafs = False
else:
  allScafs = True


if "--sep" in sys.argv:
  if getOptionValue("--sep") == "comma":
    sep = ","
  elif getOptionValue("--sep") == "white":
    sep = None
  else:
    print "\nThe only options for --sep are [comma] or [white] \n"
    sys.exit()
    
else:
  sep = None

# start output file
mainOut = open(outName, "w")
mainOut.write("scaffold,position,start,end,midpoint,sites,sitesOverMinExD")

#check analyses
analyses = []
poly = True
popPoly = False
pairWisePoly = False
polyBpp = False
popPolyBpp = False
indPoly = False
pairWisePolyBpp = False
allPairsPoly = False
if "-a" in sys.argv:
  analysesList = getOptionValue("-a").strip("\"").split(",")
  if "S" in analysesList:
    poly = True
    analyses.append("S")
    mainOut.write(",S")
  if "allPi" in analysesList:
    poly = True
    analyses.append("overallPi")
    mainOut.write(",Pi")
  if "popPi" in analysesList:
    popPoly = True
    analyses.append("popPi")
    for popName in popNames:
      mainOut.write("," + popName + "_Pi")
  if "px" in analysesList:
    popPoly = True
    analyses.append("px")
    for popName in popNames:
      mainOut.write("," + popName + "_px")
  if "indPi" in analysesList:
    indPoly = True
    analyses.append("indPi")
    for indName in indNames:
      if not includeOutGroup or indName not in vars()[outGroup]:
        mainOut.write("," + indName + "_Pi")
  if "popS" in analysesList:
    popPoly = True
    analyses.append("popS")
    for popName in popNames:
      mainOut.write("," + popName + "_S")
  if "dxy" in analysesList:
    pairWisePoly = True
    analyses.append("dxy")
    for X in range(len(popNames) - 1):
      for Y in range(X + 1,len(popNames)):
        mainOut.write("," + popNames[X] + popNames[Y] + "_dxy")
  if "ABBABABA" in analysesList:
    if "P1" in popNames and "P2" in popNames and "P3" in popNames and "O" in popNames:
      analyses.append("ABBABABA")
      mainOut.write(",ABBA,BABA,D,f,mf")
    else:
      print "\nPopulation names P1, P2, P3 and O must be present to do ABBA BABA analyses.\n"
      sys.exit()
  mainOut.write("\n")
else:
  print "\nplease specify analysis to be conducted (-a)\n"
  sys.exit()

if analyses == []:
  print "\nplease check analysis options\n"
  sys.exit()
else:
  print >> sys.stderr, "\nAnalyses to be included:\n"
  for a in analyses:
    print >> sys.stderr, a , "\n"




if "--ignoreFrequency" in sys.argv:
  iF = int(getOptionValue("--ignoreFrequency"))
else:
  iF = 0

if "--minimumExploitableData" in sys.argv:
  minExD = float(getOptionValue("--minimumExploitableData"))
  print "minimumExploitableData =", minExD
else:
  minExD = 0




# counting stat that will let keep track of how far we are
windowsTested = 0
goodWindows = 0

#create temporary  variables for nucleotide data
for name in indNames: 
  vars()["sub" + name] = []

#For the tempoarary window we need to store the positions each time to keep track of the spread of the sites 
subPos = []


#read first line and store variables
line = file.readline().rstrip()
objects = line.split(sep)
if allScafs or objects[0] in scafsOfInterest:
  subSCF = objects[0]
  subPos.append(int(objects[1]))
  for name in indNames:
    vars()["sub" + name].append(objects[names.index(name)])
else:
  subSCF = None

#read second line as the first to be evaluated by the loop
line = file.readline()
objects = line.split(sep)

windStart = 0
lastWindNA = False

while 1 == 1:
  #each time we do the loop we will be doing one window. 
  #if the line in hand is not yet too far away or on another scaffold, add the line and read another
  if allScafs or subSCF is not None:
    windowsTested += 1
  while len(objects) > 1 and objects[0] == subSCF and int(objects[1]) < windStart + windSize:
    subPos.append(int(objects[1]))
    
    for indName in indNames:
      vars()["sub" + indName].append(objects[names.index(indName)])
    
    line = file.readline()
    objects = line.split(sep)
    
  #now the line in hand is incompatible with the current window
  #if there are enough sites, we do the LD test and then slide the start along
  
  if len(subPos) >= minSites and subSCF is not None:
    
    if test or verbose:
      print "\nGood window found. Length =", len(subPos), len(vars()["sub" + indNames[0]])
    
    # add data to major outputs
    Sites = str(len(subPos))
    Scaf = (subSCF)
    Position = str(windStart + (0.5*windSize))
    Start = str(windStart)
    End = str(windStart + windSize)
    if mid(subPos):
      Midpoint = str(int(round(mid(subPos))))
    else:
      Midpoint = "NA"
    
    #if in genotype format, if diplois, split into haplos, if haploid, leave it
    if genotypes:
      for indName in indNames:
        if ploidy[indName] == 2:
          #its diploid, so split into two haplotypes
          vars()["haplo" + indName] = haplo(vars()["sub" + indName])
          vars()[indName + "A"] = vars()["haplo" + indName][::2]
          vars()[indName + "B"] = vars()["haplo" + indName][1::2]
          #if haploid, the haplotype is the same as the calls we've collected
        elif ploidy[indName] == 1:
          vars()[indName + "A"] = vars()["sub" + indName]
    # this section is for working with haplotypes separated by a | which means it must all be diploid
    elif haplotypes:
      for indName in indNames:
        vars()[indName + "A"] = []
        vars()[indName + "B"] = []
        for call in vars()["sub" + indName]:
          vars()[indName + "A"].append(call[0])
          vars()[indName + "B"].append(call[2])
    
    if test or verbose:
      print "\nHaplotypes generated. Length = ", len(vars()[indNames[1] + "A"])
    
    #create sequence objects for egglib, for all data types necessary, taking poidy into account
    #first step is to create variables for all haps and each pop which will contain a tuple for each haplotype
    allHaps = []
    for popNumber in range(len(popNames)):
      vars()[popNames[popNumber] + "Haps"] = []
      for indName in vars()[popNames[popNumber]]:
        #first, if its haploid, add only one haplotype, else add 2
        if genotypes and ploidy[indName] == 1:
          hapA = (indName + "A", "".join(vars()[indName + "A"]), popNumber)
          allHaps.append(hapA)
          vars()[popNames[popNumber] + "Haps"].append(hapA)
        else:
          hapA = (indName + "A", "".join(vars()[indName + "A"]), popNumber)
          hapB = (indName + "B", "".join(vars()[indName + "B"]), popNumber)
          allHaps.append(hapA)
          allHaps.append(hapB)
          vars()[popNames[popNumber] + "Haps"].append(hapA)
          vars()[popNames[popNumber] + "Haps"].append(hapB)
        #if we're doing any individual based heterozygosity, then make a align object for each diploid individual
        if indPoly:
          vars()[indName + "Haps"] = [hapA,hapB]
    if includeOutGroup:
      for indName in vars()[outGroup]:
        if genotypes and ploidy[indName] == 1:
          hapA = (indName + "A", "".join(vars()[indName + "A"]), 999)
          allHaps.append(hapA)
        else:
          hapA = (indName + "A", "".join(vars()[indName + "A"]), 999)
          hapB = (indName + "B", "".join(vars()[indName + "B"]), 999)
          allHaps.append(hapA)
          allHaps.append(hapB)
      #and include outgroups in populations
      for popNumber in range(len(popNames)):
        for indName in vars()[outGroup]:
          if genotypes and ploidy[indName] == 1:
            hapA = (indName + "A", "".join(vars()[indName + "A"]), 999)
            vars()[popNames[popNumber] + "Haps"].append(hapA)
          else:
            hapA = (indName + "A", "".join(vars()[indName + "A"]), 999)
            hapB = (indName + "B", "".join(vars()[indName + "B"]), 999)
            vars()[popNames[popNumber] + "Haps"].append(hapA)
            vars()[popNames[popNumber] + "Haps"].append(hapB)
    
    #now create egglib align objects for all of these sets of tuples
    # for whole set, for each pop and for pairs of pops and single inds if necessary
    allAlign = egglib.Align.create(allHaps)
    for popName in popNames:
      vars()[popName + "Align"] = egglib.Align.create(vars()[popName + "Haps"])
    if indPoly:
      for indName in indNames:
        if not includeOutGroup or indName not in vars()[outGroup]:
          vars()[indName + "Align"] = egglib.Align.create(vars()[indName + "Haps"])
        
    if pairWisePoly or pairWisePolyBpp:
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          vars()[popNames[X] + popNames[Y] + "Haps"] = []
          for hap in vars()[popNames[X] + "Haps"]:
            vars()[popNames[X] + popNames[Y] + "Haps"].append(hap)
          for hap in vars()[popNames[Y] + "Haps"]:
            vars()[popNames[X] + popNames[Y] + "Haps"].append(hap)
          vars()[popNames[X] + popNames[Y] + "Align"] = egglib.Align.create(vars()[popNames[X] + popNames[Y] + "Haps"])
    
    # if all pairwise comparisons among all individuals are being done, all pairs of two haplotypes need to be placed together in an egglib object
    if allPairsPoly:
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          for hapNumberX in range(len(vars()[popNames[X] + "Haps"])):
            for hapNumberY in range(len(vars()[popNames[Y] + "Haps"])):
              vars()[popNames[X] + popNames[Y] + "Haps" + str(hapNumberX) + str(hapNumberY)] = [vars()[popNames[X] + "Haps"][hapNumberX] , vars()[popNames[Y] + "Haps"][hapNumberY]]
              vars()[popNames[X] + popNames[Y] + "Align" + str(hapNumberX) + str(hapNumberY)] = egglib.Align.create(vars()[popNames[X] + popNames[Y] + "Haps" + str(hapNumberX) + str(hapNumberY)])
    
    if test or verbose:
      print "\negglib alignments generated:"
      print "alignment length:", allAlign.ls(), "number of sequences:", allAlign.ns()
    
    
    #depending on analyses requested, run analyses...
    if poly:
      if test or verbose:
        print "\nrunning polymorphism analyses" 
      allPoly = allAlign.polymorphism(minimumExploitableData=minExD,allowMultipleMutations=True,ignoreFrequency=iF)
    if popPoly:
      if test or verbose:
        print "\nrunning population-specific polymorphism analyses"
      for popName in popNames:
        vars()[popName + "Poly"] = vars()[popName + "Align"].polymorphism(minimumExploitableData=minExD,allowMultipleMutations=True,ignoreFrequency=iF)
    if indPoly:
      if test or verbose:
        print "\nrunning indivdual-specific polymorphism analyses"
      for indName in indNames:
        if not includeOutGroup or indName not in vars()[outGroup]:
          vars()[indName + "Poly"] = vars()[indName + "Align"].polymorphism(minimumExploitableData=minExD,allowMultipleMutations=True,ignoreFrequency=iF)
    if pairWisePoly:
      if test or verbose:
        print "\nrunning pair-wise polymorphism analyses"
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          vars()[popNames[X] + popNames[Y] + "Poly"] = vars()[popNames[X] + popNames[Y] + "Align"].polymorphism(minimumExploitableData=minExD,allowMultipleMutations=True,ignoreFrequency=iF)
    if polyBpp:
      if test or verbose:
        print "\nrunning Bio++ polymorphism analyses"
      allPolyBpp = allAlign.polymorphismBPP()
    if popPolyBpp:
      if test or verbose:
        print "\nrunning population-specific Bio++ polymorphism analyses"
      for popName in popNames:
        vars()[popName + "PolyBpp"] = vars()[popName + "Align"].polymorphismBPP()
    if pairWisePolyBpp:
      if test or verbose:
        print "\nrunning pair-wise Bio++ polymorphism analyses"
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          vars()[popNames[X] + popNames[Y] + "PolyBpp"] = vars()[popNames[X] + popNames[Y] + "Align"].polymorphismBPP()
    if allPairsPoly:
      if test or verbose:
        print "\nrunning polymorphism analyses for all pairs of individuals"
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          for hapNumberX in range(len(vars()[popNames[X] + "Haps"])):
            for hapNumberY in range(len(vars()[popNames[Y] + "Haps"])):
              vars()[popNames[X] + popNames[Y] + "Poly" + str(hapNumberX) + str(hapNumberY)] = vars()[popNames[X] + popNames[Y] + "Align" + str(hapNumberX) + str(hapNumberY)].polymorphism(minimumExploitableData=0,allowMultipleMutations=True,ignoreFrequency=iF)
    
    
    #sites passing minExD threshold
    SitesOverMinExD = str(allPoly["lseff"])
    
    #write data to main output
    mainOut.write(Scaf + "," + Position + "," + Start + "," + End + "," + Midpoint + "," + Sites + "," + SitesOverMinExD)
    
    
    if "S" in analyses:
      mainOut.write("," + str(allPoly["S"]))
      
    if "overallPi" in analyses:
      if allPoly["lseff"] >= minSites:
        mainOut.write("," + str(round(allPoly["Pi"],4)))
      else:
        mainOut.write(",NA")
    
    if "popPi" in analyses:
      for popName in popNames:
        if vars()[popName + "Poly"]["lseff"] >= minSites:
          mainOut.write("," + str(round(vars()[popName + "Poly"]["Pi"],4)))
        else:
          mainOut.write(",NA")
    
    if "px" in analyses:
      for popName in popNames:
        if vars()[popName + "Poly"]["lseff"] >= minSites:
          try:
            mainOut.write("," +  str(round(px(vars()[popName + "Align"]),4)))
          except:
            mainOut.write(",NA")
        else:
          mainOut.write(",NA")
    
    if "indPi" in analyses:
      for indName in indNames:
        if not includeOutGroup or indName not in vars()[outGroup]:
          if vars()[indName + "Poly"]["lseff"] >= minSites:
            mainOut.write("," + str(round(vars()[indName + "Poly"]["Pi"],4)))
          else:
            mainOut.write(",NA")
    
    if "popS" in analyses:
      for popName in popNames:
        if vars()[popName + "Poly"]["lseff"] >= minSites:
          mainOut.write("," + str(float(vars()[popName + "Poly"]["S"])))
        else:
          mainOut.write(",NA")
    
    if "dxy" in analyses:
      for X in range(len(popNames) - 1):
        for Y in range(X + 1,len(popNames)):
          if vars()[popNames[X] + popNames[Y] + "Poly"]["lseff"] >= minSites:
            try:
              mainOut.write("," +  str(round(dxy(vars()[popNames[X] + popNames[Y] + "Align"]),4)))
            except:
              mainOut.write(",NA")
          else:
            mainOut.write(",NA")
    
    if "ABBABABA" in analyses:
      try:
        ABstats = ABBABABA(allAlign, popNames.index("P1"), popNames.index("P2"), popNames.index("P3"), popNames.index("O"))
      except:
        ABstats = {"ABBA":"NA", "BABA":"NA", "D":"NA", "f":"NA", "mf":"NA"}
      mainOut.write("," + str(ABstats["ABBA"]))
      mainOut.write("," + str(ABstats["BABA"]))
      mainOut.write("," + str(ABstats["D"]))
      mainOut.write("," + str(ABstats["f"]))
      mainOut.write("," + str(ABstats["mf"]))
    
    
    mainOut.write("\n")
    
    
    goodWindows += 1
    
    if test:
      break
    
    
    if stopAt:
      if stopVal == goodWindows:
        break
    
    lastWindNA = False
    
    windStart += slide
    i = len(subPos)
    for x in subPos:
      if x >= windStart:
        i = subPos.index(x)
        break
    subPos = subPos[i:]
    
    for name in indNames:
      vars()["sub" + name] = vars()["sub" + name][i:]
  
  #otherwise, if the last window as not NA, we will make an NA window and then we'll slide along (or reset if we're onto a new scaf
  else:
    if subSCF is not None and lastWindNA == False:
      Sites = str(len(subPos))
      SitesOverMinExD = "NA"
      Scaf = (subSCF)
      Position = str(windStart + (0.5*windSize))
      Start = str(windStart)
      End = str(windStart + windSize)
      if mid(subPos):
        Midpoint = str(int(round(mid(subPos))))
      else:
        Midpoint = "NA"
      
      mainOut.write(Scaf + "," + Position + "," + Start + "," + End + "," + Midpoint + "," + Sites + "," + SitesOverMinExD)
      
      #Fill in NAs for all requested data
      if "overallPi" in analyses:
        mainOut.write(",NA")
      
      if "popPi" in analyses:
        for popName in popNames:
          mainOut.write(",NA")
      
      if "px" in analyses:
        for popName in popNames:
          mainOut.write(",NA")
      
      if "indPi" in analyses:
        for indName in indNames:
          if not includeOutGroup or indName not in vars()[outGroup]:
            mainOut.write(",NA")
      
      if "popS" in analyses:
        for popName in popNames:
          mainOut.write(",NA")
      
      if "dxy" in analyses:
        for X in range(len(popNames) - 1):
          for Y in range(X + 1,len(popNames)):
            mainOut.write(",NA")
      
      if "ABBABABA" in analyses:
        mainOut.write(",NA,NA,NA,NA,NA")
      
      if "S" in analyses:
        mainOut.write(",NA")
      
      #and end the line
      mainOut.write("\n")
      
      
      #and record the winow as an NA
      lastWindNA = True
    
    #if the line in hand is on the same scaf, we simply slide the start along one slide
    if len(objects) > 1 and objects[0] == subSCF:
      i = len(subPos)
      windStart += slide
      for x in subPos:
        if x >= windStart:
          i = subPos.index(x)
          break
      subPos = subPos[i:]
      
      for name in indNames:
        vars()["sub" + name] = vars()["sub" + name][i:]
    #otherwise its a new scaf, so we reset the subwindow and subScaf and read the next line
    else:
      windStart = 0
      
      if len(objects) > 1:
        if allScafs or objects[0] in scafsOfInterest:
          subSCF = objects[0]
          subPos = [int(objects[1])]
      
          for name in indNames:
            vars()["sub" + name] = [objects[names.index(name)]]
        else:
          subSCF = None
      
        line = file.readline().rstrip()
        objects = line.split(sep)

      else:
        break
  if windowsTested > 0 and windowsTested % report == 0:
    print windowsTested, "windows done ..."


file.close()
mainOut.close()


print "\n" + str(windowsTested) + " windows were tested."
print "\n" + str(goodWindows) + " windows were good.\n"

print "\nDone."
