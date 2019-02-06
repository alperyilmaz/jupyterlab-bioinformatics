import itertools as it
import sys
import random as rd
import networkx as nx
from collections import defaultdict
from collections import Counter


def redundantSimple(G,seeds,tech,techargs):
	Gs=tech(G,seeds,techargs,1)
	numNodes=Gs.number_of_nodes()
	numEdges=Gs.number_of_edges()
	redundantSeeds=[]
	for i in range(len(seeds)):
		Gs1=tech(G,seeds[:i]+seeds[i+1:],techargs,1)
		if Gs1.number_of_nodes()==numNodes and Gs1.number_of_edges()==numEdges:
			redundantSeeds.append(seeds[i])
	remainingSeeds=[x for x in seeds if x not in redundantSeeds]
	foundSolution=0
	solution=[]
	for i123 in range(len(redundantSeeds)):
		for combin in it.combinations(redundantSeeds,i123):
			Gs1=tech(G,remainingSeeds+list(combin),techargs,1)
			if Gs1.number_of_nodes()==numNodes and Gs1.number_of_edges()==numEdges:
				foundSolution=1
				solution.append(remainingSeeds+list(combin))
	#		print solution
		if foundSolution==1:
			return solution
	return [seeds]
