import networkx as nx

def globalmetrics(G):
	result={}
	result['num_nodes']=G.number_of_nodes()
	result['num_edges']=G.number_of_edges()
	try:
		result['Average Degree']=float(sum(G.degree().values()))/float(G.number_of_nodes())
	except:
		result['Average Degree']=0
	try:
		result['clustering']=nx.average_clustering(G)
	except:
		result['clustering']=0
	try:
		result['assortivitity']= assortativity_nx(G) #nx.assortativity.degree_assortativity_coefficient(G)
	except:
		result['assortivitity']=0
	result['components']=nx.number_connected_components(G)
	return result

def assortativity_nx(graph):
    # implementation taken and optimised from: https://github.com/dougvk/CS462/blob/master/hw1/measure_graph.py
    # comment from original document detailing where this came from: igraph doesn't have calculation -- found one online @ snipplr.com
    degrees = graph.degree()
    degrees_sq = {deg:degrees[deg]**2 for deg in degrees}

    m = float(graph.number_of_edges())
#    num1, num2, den1 = 0, 0, 0
#    for source, target in graph.get_edgelist():
#        num1 += degrees[source] * degrees[target]
#        num2 += degrees[source] + degrees[target]
#        den1 += degrees_sq[source] + degrees_sq[target]
    edgelist=graph.edges()
    num1=sum([degrees[source]*degrees[target] for source, target in edgelist])
    num2=sum([degrees[source]+degrees[target] for source, target in edgelist])
    den1=sum([degrees_sq[source]+degrees_sq[target] for source, target in edgelist])

    num1 /= m
    den1 /= 2*m
    num2 = (num2 / (2*m)) ** 2

    return (num1 - num2) / (den1 - num2)
