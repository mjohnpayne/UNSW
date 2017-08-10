#!/usr/bin/python 

import sys
import re 
import logging
# logging.basicConfig(level=logging.DEBUG)


def printLineOutput(seqId, location, score, seqName): 
 	print seqId + "\t" + location + "\t" + score + "\t" + seqName; 


def getFinalPredictions(filenameIn, cyto, cytoMembrane, peri, outerMembrane, extracellular, unknown): 

	fh = open(filenameIn, 'r') 



	seqId = ""
	flag_isNextFinal = False; 

	print("SeqId\tLocalisation\tScore\tSeqenceName") 

	for line in fh: 
		line = line.strip()

		if re.match("^SeqID:", line): 
			arr = re.split("\|", line)

			seqId = arr[1] 
			seqName = arr[4]

			logging.debug("line: " + line) 
			logging.debug("seqId: " + arr[1] ) 

			flag_isNextFinal = False; 

		if re.match("^final", line, flags=re.IGNORECASE): 
			flag_isNextFinal = True; 

		elif flag_isNextFinal == True: 
			arr = re.split("[\s\t]+", line) 

			if cyto != "" and re.match("^" + cyto + "[\s\t]+", line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], arr[len(arr)-1], seqName) 
			


			if cytoMembrane != "" and re.match("^" + cytoMembrane + "[\s\t]+", line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], arr[len(arr)-1], seqName)


			if peri != "" and re.match("^" + peri + "[\s\t]+", line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], arr[len(arr)-1], seqName)

			if outerMembrane != "" and re.match("^" + outerMembrane + "[\s\t]+", line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], arr[len(arr)-1], seqName)

			if extracellular != "" and re.match("^" + extracellular + "[\s\t]+", line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], arr[len(arr)-1], seqName)

			if unknown != "" and re.match("^" + unknown, line, flags=re.IGNORECASE): 
				printLineOutput(seqId, arr[0], "-", seqName)

			flag_isNextFinal = False

		
	 


def getOptions(args): 
	filename = "" 
	cyto = ""
	cytoMembrane = "" 
	peri = ""
	outerMembrane = ""
	extracellular = ""
	unknown = "" 

	for i in range(1, len(args)): 
		args[i] = re.sub('[\s\t]+', '', args[i], flags=re.IGNORECASE)

		if re.match("^[^-]", args[i], flags=re.IGNORECASE): 
			filename = args[i] 
			logging.debug('filename: ' + filename) 

		if re.match("^-cytoplasmic$", args[i], flags=re.IGNORECASE): 
			cyto = re.sub("-", "", args[i]) 
			logging.debug('cyto: ' + cyto) 


		if re.match("^-cytoplasmicmembrane$", args[i], flags=re.IGNORECASE):
			cytoMembrane =  re.sub("-", "", args[i]) 
			logging.debug('cytoMembrane: ' + cytoMembrane) 

		if re.match("^-periplasmic$", args[i], flags=re.IGNORECASE):
			peri =  re.sub("-", "", args[i]) 
			logging.debug('peri: ' + peri) 

		if re.match("^-outermembrane$", args[i], flags=re.IGNORECASE):
			outerMembrane =  re.sub("-", "", args[i]) 
			logging.debug('outerMembrane: ' + outerMembrane) 

		if re.match("^-extracellular$", args[i], flags=re.IGNORECASE):
			extracellular =  re.sub("-", "", args[i]) 
			logging.debug('extracellular: ' + extracellular) 

		if re.match("^-unknown$", args[i], flags=re.IGNORECASE):
			unknown =  re.sub("-", "", args[i]) 
			logging.debug('unknown: ' + unknown) 
		

	if filename == "": 
		print "Error: input filename not supplied"
		printErrorAndExit()
	


	return (filename, cyto, cytoMembrane, peri, outerMembrane, extracellular, unknown); 




def printErrorAndExit(): 
	usage = "python script.py <localisation.txt> [-cytoplasmic -cytoplasmicmembrane -periplasmic -outermembrane -extracellular -unknown]\n\r" 
	print "Usage: " + usage; 
	exit() 



if __name__ == '__main__':

	if len(sys.argv) < 3 or len(sys.argv) > 8: 
		print "Error: incorrect number of inputs\n\r"
		printErrorAndExit(); 

	(filename, cytoplasmic, cytoMembrane, peri, outerMembrane, extracellular, unknown) = getOptions(sys.argv) 


 	getFinalPredictions(sys.argv[1], cytoplasmic, cytoMembrane, peri, outerMembrane, extracellular, unknown) 












