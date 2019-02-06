import sys
import networkx as nx
import itertools as it
precomputedPaths=[]
from collections import OrderedDict




## Snowball Sampling Code


def __snowballsample2hop_igraph__(G,seeds): #tested
    seeds1=seeds
    seeds2=[]
    for x in seeds1:
        seeds2+=G.neighbors(x)
    seeds2=list(set(seeds2))
    for x in seeds2:
        seeds1+=G.neighbors(x)
    seeds2+=seeds1
    return G.subgraph(list(set(seeds2)))




def __snowballsample__(G,node_list,number_of_hops=2): #testing code
    nodes2add=[list(node_list)]
    for hop in range(number_of_hops):
        temp=[]
        for node in nodes2add[hop]:
            temp+=G[node]
        nodes2add.append(temp)
    H=G.subgraph(sum(nodes2add,[]))
    H.to_undirected()
    return H


def snowballsampling(G,seeds,hops=2):
    if str(type(G))=="<class 'networkx.classes.graph.Graph'>":
        return __snowballsample__(G,seeds,hops)
    elif str(type(G))=="<class 'igraph.Graph'>":
        if hops==2:
            return __snowballsample2hop_igraph__(G,seeds)
        else:
            print "Method for snowball sampling igraph networks has not been implemented"
            assert(0==1)


def __pred_maker__(G,source): #testing code
    level=0                  # the current level
    nextlevel=[source]       # list of nodes to check at next level
    seen={source:level}      # level (number of hops) when seen in BFS
    pred={source:[]}         # predecessor dictionary
    while nextlevel:
        level=level+1
        thislevel=nextlevel
        nextlevel=[]
        for v in thislevel:
            for w in G[v]:
                if w not in seen:
                    pred[w]=[v]
                    seen[w]=level
                    nextlevel.append(w)
                elif (seen[w]==level):# add v to predecessor list if it
                    pred[w].append(v) # is at the correct level
    return pred

def __pred_maker_two__(G,source): #testing code
    level=0                  # the current level
    nextlevel=[source]       # list of nodes to check at next level
    seen={source:level}    # level (number of hops) when seen in BFS
    pred={x:[-1] for x in G.nodes()}
    pred[source]=[]	    # predecessor dictionary
    while nextlevel:
        level=level+1
        thislevel=nextlevel
        nextlevel=[]
        for v in thislevel:
            for w in G[v]:
                if w not in seen:
                    pred[w]=[v]
                    seen[w]=level
                    nextlevel.append(w)
                elif (seen[w]==level):# add v to predecessor list if it
                    pred[w].append(v) # is at the correct level
    return pred


def __pred_smaller__(G,source):
	print source
	pred1=__pred_maker_two__(G,source)
	result1=[tuple(pred1[x]) for x in range(len(pred1))]
	result=tuple(result1)
	return result



def __G_whole_pred__(prot):
	return __pred_maker__(rb.G_whole,prot)





def shortest_paths_sampling(graph,seedlist1):
    if str(type(graph))=="<class 'networkx.classes.graph.Graph'>":
        result=__shortest_paths_sampling_networkx__(graph,seedlist1)
    elif str(type(graph))=="<class 'igraph.Graph'>":
        result=__shortest_paths_sampling_igraph__(graph,seedlist1)
    else:
        import sys
        sys.exit()
    return result

def __make_precomputed_paths__():
	import gc
	graph=rb.G_whole
	result={}
	count=0
	for protein in graph.nodes():
		result[rb.no[protein]]={}
		count+=1
		print count
		pred=__pred_maker__(graph,protein)
		result[protein]=pred
#		for item in pred:
#			if rb.no[item]> rb.no[protein]:
#				result[(rb.no[protein],rb.no[item])]=tuple(set(__shortest_path_recursion_edgePairs__(pred,item)))
		if count%100==0:
			gc.collect()
#			return result
	return result


def test_shortestpath(graph,seedlist):
	result={}
	G=nx.Graph()
	for i1 in range(len(seedlist)-1):
		protein=seedlist[i1]
		pred=preds[protein]
		test=[]
		for i2 in range(i1+1,len(seedlist)):
			test+=__shortest_path_recursion_edgePairs__(pred,seedlist[i2])
		G.add_edges_from(test)
	return G







def __shortest_path_recursion_edgePairs__(pred,endpoint):
    if endpoint not in pred:
        print 'endpoint not in pred'
        return []
    if pred[endpoint]==[]:
        return []
    listofEdges=[]
    pos=endpoint
    while len(pred[pos])==1:
	listofEdges+=[(rb.no[pos],rb.no[pred[pos][0]])]
	pos=pred[pos][0]
    if pred[pos]==[]:
	return listofEdges
    else:
        for item in pred[endpoint]:
            listofEdges+=__shortest_path_recursion_edgePairs__(pred,item)
    return listofEdges


def __shortest_path_recursion_edgePairs_production__(pred,endpoint):
    if endpoint not in pred:
        print 'endpoint not in pred'
        return []
    if pred[endpoint]==[]:
        return []
    listofEdges=[]
    pos=endpoint
    while len(pred[pos])==1:
       listofEdges+=[(pos,pred[pos][0])]
       pos=pred[pos][0]
    if pred[pos]==[]:
        return listofEdges
    else:
        for item in pred[pos]:
            listofEdges+=[(pos,item)]
            listofEdges+=__shortest_path_recursion_edgePairs_production__(pred,item)
    return listofEdges




def __get_precomputed_paths__():
	return []

def __shortestpathPrecomputed__(seedlist1):
	if precomputedPaths==[]:
		precomputedPaths=__get_precomputed_paths__()
	G=nx.Graph()
	G.add_edges_from(precomputedPaths[frozenset(x)] for x in it.combination(seedlist1))
	return G



def __shortest_paths_sampling_igraph__(graph,seeds):
    G=nx.Graph()
    te=[]
    for seed in seeds:
    	te+=[x for x in graph.get_all_shortest_paths(seed) if len(x)>1 if x[-1] in seeds]

    G.add_nodes_from(seeds)
    G.add_edges_from((x[i],x[i+1]) for i in range(len(x)-2) for x in te)

    return G








def __shortest_paths_sampling_networkx_paths__(G,seedlist): #Tested
	paths=[]
	for i in range(len(seedlist)):
		seed=seedlist[i]
		pred = __pred_maker__(G,seed)
		for j  in range(i,len(seedlist)):
			if seedlist[j] in pred:
#				print (seedlist[i],seedlist[j])
				paths+=__shortest_path_recursion_edgePairs_production__(pred,seedlist[j])

#		else:
#		print str(seedlist[j]) +' unreachable from '+str(seedlist[i])
	return paths

def __shortest_paths_sampling_networkx__(G,seedlist): #Tested
    paths=__shortest_paths_sampling_networkx_paths__(G,seedlist)
    G=nx.Graph()
    G.add_edges_from(paths)
    G.to_undirected()
    return G

def __shortest_paths_sampling_pathsInfo__(G,seedlist): #Tested
    paths=__shortest_paths_sampling_networkx_paths__(G,seedlist)
    G=nx.Graph()
    G.add_edges_from(paths)
    G.to_undirected()
    return (G,paths)


def __shortest_paths_sampling_networkx_old__(G,seedlist):
    paths=[]
    #max1={x:max([nx.shortest_path_length(G,x,y) for y in seedlist]) for x in seedlist}
    for i in range(len(seedlist)):
        seed=seedlist[i]

        pred = __pred_maker_old__(G,seed)#,max1[seed])
        paths+=[(t[i],t[i+1]) for j  in range(i,len(seedlist)) for t in __shortest_path_recursion__(pred,seedlist[j]) for i in range(len(t)-1)]

    G=nx.Graph()
    G.add_edges_from(paths)
    G.to_undirected()
    return G


def __shortest_path_recursion__(pred,endpoint):
    if pred[endpoint]==[]:
        return [[endpoint]]
    listofpaths=[]
    for item in pred[endpoint]:
        listofpaths+=__shortest_path_recursion__(pred,item)

    for path in listofpaths:
        path.append(endpoint)

    return listofpaths



def all_paths_shorter_than_k(graph,seeds,k=4):
    if str(type(graph))=="<class 'networkx.classes.graph.Graph'>":
        result=__all_paths_networkx__(graph,seeds,k)
    elif str(type(rb.G_ig))=="<class 'igraph.Graph'>":
        result=__all_paths_igraph__(graph,seeds,k)
    else:
    	print 'Not a valid graph file'
        import sys
        sys.exit()
    return result


def __paths_less_than_k_networkx_recursion_new__(graph,node,cur_path,level,k,path_dict,path_dict_path,nSeed):
    if node not in path_dict:
        path_dict[node]=[nSeed]
    elif nSeed not in path_dict[node]:
        path_dict[node].append(nSeed)
    if node not in path_dict_path:
        path_dict_path[node]=[cur_path]
    else:
        path_dict_path[node].append(cur_path)

    path=[cur_path]
    if level!=k:
        for n in graph[node]:
            __paths_less_than_k_networkx_recursion_new__(graph,n,cur_path+[n],level+1,k,path_dict,path_dict_path,nSeed)

    return



def __all_paths_networkx__(graph,seeds,k):

    if k==4:
    	path=[]
    	endpoint2seed={}
    	for seed in seeds:
    		for n1 in graph[seed]:
    			path.append([seed,n1])
    			if n1 in endpoint2seed:
    				endpoint2seed[n1].append(seed)
    			else:
     				endpoint2seed[n1]=[seed]
    			for n2 in graph[n1]:
    				path.append([seed,n1,n2])
    				if n2 in endpoint2seed:
    					endpoint2seed[n2].append(seed)
    				else:
     					endpoint2seed[n2]=[seed]

     	mask1={x for x in endpoint2seed if len(set(endpoint2seed[x]))>1}

     	G=nx.Graph()
     	G.add_edges_from((x[i],x[i+1]) for x in path if x[-1] in mask1 for i in range(len(x)-1))
     	return G
    return __all_paths_networkx_helper__(graph,seeds,k)


def __all_paths_networkx_helper__(graph,seeds,k):
    result=[]
    paths={}
    path_dict={}
    path_dict_path={}

    if k==1:
 #       print 'used k=1 speedup ' +str(seeds)

	G=nx.Graph()
	G.add_edges_from(x for x in it.ifilter(lambda x: graph.has_edge(x[0],x[1]),it.combinations(seeds,2)))
	return G


    if k%2==0:
        for seed in seeds:
            cur_path=[seed]
            __paths_less_than_k_networkx_recursion_new__(graph,seed,cur_path,0,k/2,path_dict,path_dict_path,seed)

        count=0

        G=nx.Graph()
        G.add_edges_from(iter(generator_test(path_dict,path_dict_path)))
        return G
    else:
  #      print 'i used the new odd code'

        return __all_paths_odd_seeds__(graph,seeds,k)


        return __all_paths_networkx_slow__(graph,seeds,k)
#    for endpoint in path_dict:
#        t=[]
#        if len(path_dict[endpoint])>1:
#            for cur_path in path_dict_path[endpoint]:
#                G.add_edges_from(pairwise(iter(cur_path)))
                #G.add_edges_from([(cur_path[i],cur_path[i+1]) for i in range(len(cur_path)-1)])
        #G.add_edges_from(t)
#    return G


def __all_paths_odd_seeds__(graph,seeds,k):
    result=[]
    paths={}
    path_dict={}
    path_dict_path={}
    path_dict_odd={}
    path_dict_path_odd={}
    for seed in seeds:
        cur_path=[seed]
        __paths_less_than_k_networkx_recursion_odd__(graph,seed,cur_path,0,(1+k)/2,path_dict,path_dict_path,seed,path_dict_odd,path_dict_path_odd)
    count=0
    G=nx.Graph()
    G.add_edges_from(iter(generator_test(path_dict,path_dict_path)))
#    print "Im before the loop"
#    print path_dict_odd
#    return (path_dict,path_dict_path,path_dict_odd,path_dict_path_odd)
#    return G
#    print 'Im adding edges now'
    for x in path_dict_odd:
        if x in path_dict:
             if len(set(path_dict_odd[x]+path_dict[x]))==1:
                continue
             if len(path_dict[x])>1:
#                if len(path_dict[x])==1:
#                   G.add_edges_from([(y[i],y[i+1]) for y in path_dict_path_odd[x]+path_dict_path[x] for i in range(len(y)-1)])
#                else:
                G.add_edges_from((y[i],y[i+1]) for y in path_dict_path_odd[x] for i in range(len(y)-1))
             else:
                exc=path_dict[x][0]
                G.add_edges_from((y[i],y[i+1]) for y in path_dict_path_odd[x] if y[0]!=exc for i in range(len(y)-1))
                G.add_edges_from((y[i],y[i+1]) for y in path_dict_path[x] for i in range(len(y)-1))


    return G



def __paths_less_than_k_networkx_recursion_odd__(graph,node,cur_path,level,k,path_dict,path_dict_path,nSeed,path_dict_odd,path_dict_path_odd):
#    print 'odd recursion'

    if node not in path_dict:
        path_dict[node]=[nSeed]
    elif nSeed not in path_dict[node]:
        path_dict[node].append(nSeed)

    if node not in path_dict_path:
        path_dict_path[node]=[cur_path]
    else:
        path_dict_path[node].append(cur_path)

    path=[cur_path]
    if level!=k-1:
        for n in graph[node]:
            __paths_less_than_k_networkx_recursion_odd__(graph,n,cur_path+[n],level+1,k,path_dict,path_dict_path,nSeed,path_dict_odd,path_dict_path_odd)
    elif level==k-1: # and k%2==1:
	n=node
        cur_path1=cur_path[:]
        for node in graph[n]:
           if node==cur_path1[0]:
		continue
           cur_path=cur_path1+[node]
           if node not in path_dict_odd:
              path_dict_odd[node]=[nSeed]
           elif nSeed not in path_dict_odd[node]:
              path_dict_odd[node].append(nSeed)

           if node not in path_dict_path_odd:
              path_dict_path_odd[node]=[cur_path]
           else:
              path_dict_path_odd[node].append(cur_path)
    return











def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = it.tee(iterable)
    next(b, None)
    return it.izip(a, b)



def generator_test(path_dict,path_dict_path):
    for endpoint in path_dict:
        t=[]
        if len(path_dict[endpoint])>1:
            for cur_path in path_dict_path[endpoint]:
                for i in range(len(cur_path)-1):
                    yield (cur_path[i],cur_path[i+1])
    return







def __all_paths_networkx_slow__(graph,seeds,k):
    paths=[]
    for seed in seeds:
        cur_path=[seed]
        paths+=__paths_less_than_k_networkx_recursion_slow__(graph,seed,cur_path,0,k,seeds)
   # return paths
    G=nx.Graph()
    for path in paths:
        if len(path)<2:
           continue
        assert(path[0] in seeds)
        assert(path[-1] in seeds)
        temp=[(path[i],path[i+1]) for i in range(len(path)-1)]
        G.add_edges_from(temp)
    return G







def __paths_less_than_k_networkx_recursion_slow__(graph,node,cur_path,level,k,seeds):
    path=[]
    if cur_path[-1] in seeds:
        if cur_path[-1]==cur_path[0] and len(cur_path)>1:
            return path
        path+=[cur_path[:]]
    if level!=k:
        for n in graph[node]:
            temp=__paths_less_than_k_networkx_recursion_slow__(graph,n,cur_path+[n],level+1,k,seeds)
            path+=temp[:]
    return path










## Simple paths less than or equal to k
from collections import deque
def simplePathlessthank(G,seeds,k):

	G1=nx.Graph()
	#We are going to make this as simple as possible
	seeds.sort(key=lambda x : len(G[x]))
	Q=deque([[x,] for x in seeds[:-1]])
	while Q:
		cur=Q.pop()
		len1=len(cur)
		if len1==k:
			l1=[n for n in G[cur[-1]] if n in seeds and n!=cur[0]]
			if len(l1)>0:
				G1.add_edges_from(((cur[i],cur[i+1]) for i in range(len1-1)))
				G1.add_edges_from(((cur[-1],n) for n in l1))
		else:

			for n in G[cur[-1]]:  #if x not in cur]:
				cur1=cur+[n]
#				print cur1
				if n in cur:
					continue
				if n in seeds and n!=cur1[0]:
					G1.add_edges_from(((cur1[i],cur1[i+1]) for i in range(len1)))
					continue
				if (len(cur1)<k+1): # and n not in seeds):
					Q.append(cur1)
	return G1

def simplePathIgraph(G,seeds,k,n1=0):
	if n1==1:
		return simplePathlessthank(G,seeds,k)
	return rb.networkx_to_igraph(simplePathlessthank(G,seeds,k))
