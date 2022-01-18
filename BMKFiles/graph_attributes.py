import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import csv
from collections import Counter

#Use this script to extract graph attributes for your parsed models list (returned from parse_convert.py)

def import_graphs(file):
	line = file.readline()
	roots = []
	master_edges = []
	while line:
		if line[0].isdigit(): #found model

			line = file.readline() #root line
			for s in line.split(): 
				if s.isdigit():
					roots.append(s)
			line = file.readline() #first digit
			num_edges = 0
			edges = []
			while line[0].isdigit():
				num_edges = num_edges +1
				for s in line.split(): 
					if s.isdigit():
						edges.append(s)
				line = file.readline()
			it = iter(edges)
			edges = list(zip(it,it))
			#print(edges)
			master_edges.append(edges)
		line = file.readline()
	return roots, master_edges;
	

def degree_of_root(roots,master_edges):
	root_degrees = []
	for i in range(len(roots)): 
		root_count = 0
		for j in master_edges[i]:
			if j[0] == roots[i] or  j[1] == roots[i]:
				root_count = root_count + 1
				
		root_degrees.append(root_count)		
	return(root_degrees)
	
	
def num_edges(master_edges):
	n_edges = []
	for i in range(len(master_edges)): 
		
				
		n_edges.append(len(master_edges[i]))		
	return(n_edges)
	
def mean_shortest_path_length(master_edges,roots):
	msp = []
	
	for e in range(0,len(master_edges)):
		G=nx.Graph(master_edges[e])
		lengths = nx.single_source_shortest_path_length(G,roots[e]) #returns a dictionary
		print(lengths)
		counter = 0
		for key in lengths:
			counter = counter + lengths[key]
		msp.append(counter/len(G))
	return(msp)	
	
def mean_clustering_coeff(master_edges):
	mcc = []
	for e in range(0,len(master_edges)):
		G=nx.Graph(master_edges[e])
		print(nx.average_clustering(G))
		mcc.append(nx.average_clustering(G))
	return(mcc)	
		
def find_cliques(master_edges,roots):
	cliq = []
	for e in range(0,len(master_edges)):
		G=nx.Graph(master_edges[e])
		print(nx.node_clique_number(G,roots[e]))
		
	return(cliq)			
	
def find_cycle_length_with_root(master_edges,roots):
	rootin3_4 = []
	for e in range(0,len(master_edges)):
		G=nx.Graph(master_edges[e])
		cycle_basis = nx.cycle_basis(G)
		#print(cycle_basis)
		
		if len(cycle_basis) > 0:
			found_cycle = 0
			for cycle in cycle_basis:
				if roots[e] in cycle:
					if len(cycle) == 3 or len(cycle) == 4:
						found_cycle = 1
						break
			if found_cycle == 1:
				rootin3_4.append(1)
			else:
				rootin3_4.append(0)
		else:
			rootin3_4.append(0)
	return(rootin3_4)	

def contains_cycle(master_edges):
	cycles = []
	for e in range(0,len(master_edges)):
		G=nx.Graph(master_edges[e])
		cycle_basis = nx.cycle_basis(G)
		#print(cycle_basis)
		
		if len(cycle_basis) > 0:
			cycles.append(1)
		else:
			cycles.append(0)
	return(cycles)
	
class MarkovGraph:
	def __init__(self, N, edgelist, root):
		self.N = N
		self.edges = list(edgelist)
		self.root = root
		self.adjm = np.zeros((N,N))
	def build_adj(self):
		for edge in self.edges:
			self.adjm[edge[0]][edge[1]] = 1;
			self.adjm[edge[1]][edge[0]] = 1;
		#print(self.adjm)
	def find_max_degree(self):
		max_degree = 0
		for i in range(self.N):
			degree = sum(self.adjm[i])
			if degree > max_degree:
				max_degree = degree
		return max_degree
	def print_cycles(self):	
		G=nx.Graph(self.edges)
		return(list(nx.cycle_basis(G)))
	def print_topology(self,fparsed):
		fparsed.write(f'Root:\t{self.root}\n')
		for edge in self.edges:
			fparsed.write('\t'.join(str(s) for s in edge)+'\n')
		fparsed.write('\n')
		
		
def main():
	file = open(sys.argv[1], "r") # your StateNparsedDT"X"CLT"Y".txt file returned from parse_convert.py
	roots, master_edges = import_graphs(file)
	
	# msp = mean_shortest_path_length(master_edges,roots)
	# mcc = mean_clustering_coeff(master_edges)
	#cliq = find_cliques(master_edges,roots)
	edge_numbers = num_edges(master_edges)
	root_degrees = degree_of_root(roots,master_edges)
	cycles = contains_cycle(master_edges)
	cycleswroot = find_cycle_length_with_root(master_edges,roots)
	
	with open('ModelAttributes.csv', mode='w',newline = '') as csv_file:
		fieldnames = ['ID', 'N_Edges', 'Root Degree','Contains Cycle','Root in 3 or 4 Cycle']
		writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

		writer.writeheader()
		for i in range(len(roots)):
			writer.writerow({'ID': i+1, 'N_Edges': edge_numbers[i], 'Root Degree': root_degrees[i], 'Contains Cycle': cycles[i], 'Root in 3 or 4 Cycle': cycleswroot[i]})
		
	
if __name__ == "__main__":
	main()
