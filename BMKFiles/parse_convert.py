# File: parse_and_convert_models.py



##This script takes in the Model"N"BMK.txt files and returns a parsed model list (StateNparsedDT"X"CLT"Y".txt with specified degree (X) and cycle length (Y) restrictions.


import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


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
		
		
def convert_and_parse(file, N,degree_threshold,cycle_length_threshold):
	fparsed = open("State" + str(N)+ "parsed" + "DT" + str(degree_threshold) + "CLT" + str(cycle_length_threshold)+ ".txt", "w+")
	#num_edges = open("State" + str(N) + "edges.txt","w+")
	line = file.readline()
	counter_include = 0
	counter = 0
	exclude_indices = []
	while line:
		counter = counter + 1 
		if counter % 10000 == 0:
			print("Processing Model:\t" + str(counter))
		if N < 10:
			start = 4
			pos = 2 + N*2 + 1+1
		else:
			start  = 5; 
			pos = 3 + N*2 + 1+1; 
		parsed_line = line[start:pos]	
		#print(parsed_line)
		root = parsed_line.split().index("1")
		parsed_line1 = line[pos:]
		edges = parsed_line1.split()
		edges = list(map(int,edges))
		it = iter(edges)
		edges = zip(it,it)
		Graph = MarkovGraph(N,edges,root)
		#print(counter)
		#print('Root:\t'+ str(root))
		Graph.build_adj()
		max_degree = Graph.find_max_degree()
		#print(max_degree)
		if max_degree <= degree_threshold:
			cycles = Graph.print_cycles()
			#print(cycles)
			if (len(cycles) > 0):
				max_cycle_length = max(len(cycle) for cycle in cycles)
			else:	
				max_cycle_length = 0
			if(max_cycle_length <= cycle_length_threshold):
				#topology passes
				counter_include = counter_include+1
				fparsed.write(f'{counter}\n')
				Graph.print_topology(fparsed)
			else: #model passes degree reqs but not cycle length reqs
				exclude_indices.append(counter)
		else: #model already does not pass due to degree requirements
			exclude_indices.append(counter)
		line = file.readline()
	fparsed.write(f'Count Included Graphs:\t{counter_include}\n')
	#print("Count Included Graphs:\t" + str(counter_include) + "\n")
	fparsed.write("not viable indices\n")
	fparsed.write(str(exclude_indices))	
	fparsed.close()
	
		
        

def main():
	file = open(sys.argv[1], "r") #corresponding BMK file for desired number of states to parse
	convert_and_parse(file,5,4,4) #number of states, degree restriction, cycle length restriction
	file.close()
	
if __name__ == "__main__":
	main()

