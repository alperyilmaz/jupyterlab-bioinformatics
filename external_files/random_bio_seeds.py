import math
import time
import networkmetricsets as ms
import random as rd
import sampling as sp
import inspect
import sys
from multiprocessing import Pool
try:
	import igraph as ig
except:
	print 'i could not import igraph'
import networkx as nx


###### Sampling Techniques Wrappers #######

def samplePathsShorterThan(G,seeds,k=4,networkx=0,markseeds=0): #tested by individual functions
	if networkx==0:
		G1=networkx_to_igraph(sp.all_paths_shorter_than_k(G,seeds,k))
		return G1
	else:
		if markseeds==0:
			return sp.all_paths_shorter_than_k(G,seeds,k)
		else:
			return __addseeds__(sp.all_paths_shorter_than_k(G,seeds,k),seeds)


def sampleShortestPaths(G,seeds,extraArgs=[],networkx=0,markseeds=0): #tested by individual functions
	if networkx==0:
		G1=networkx_to_igraph(sp.shortest_paths_sampling(G,seeds))
		return G1
	else:
		if markseeds==0:
			return sp.shortest_paths_sampling(G,seeds)
		else:
			return __addseeds__(sp.shortest_paths_sampling(G,seeds),seeds)

def sampleSnowball(G,seeds,hops=[],networkx=0):
	if networkx==1:
		return sp.snowballsampling(G,seeds,hops)
	if G.__class__==nx.classes.graph.Graph:
		if networkx==0:
			return networkx_to_igraph(sp.snowballsampling(G,seeds,hops))
		else:
			return sp.snowballsampling(G,seeds,hops)

	if hops==2:
		seeds1=seeds[:]
		seeds2=[]
		for x in seeds1:
			seeds2+=G.neighbors(x)
		seeds2=list(set(seeds2))
		for x in seeds2:
			seeds1+=G.neighbors(x)
		seeds2+=seeds1
		seeds2=list(set(seeds2))
		return G.subgraph(seeds2)
	elif hops==1:
		seeds1=seeds[:]
		seeds2=[]
		for x in seeds1:
			seeds2+=G.neighbors(x)
		seeds2=list(set(seeds2+seeds1))
		return G.subgraph(seeds2)
	else:
		return networkx_to_igraph(sp.snowballsampling(G,seeds,hops))









###### Randomisations #######


def __make_num_dict__(G): #tested
    deg_list={x:len(G[x]) for x in G}
    num_deg_dict={}
    for item in deg_list:
        if deg_list[item] in num_deg_dict:
            num_deg_dict[deg_list[item]].append(item)
        else:
            num_deg_dict[deg_list[item]]=[item]
    return num_deg_dict


def __make_bins__(data_dict,min_num=5): #tested
    resulting_bins={}
    kys=data_dict.keys()
    kys.sort()
    temp=[]
    key_set=[]
    for ky in kys:
        temp+=data_dict[ky]
        key_set.append(ky)
        if len(temp)>min_num-1:
            for ky1 in key_set:
                resulting_bins[ky1]=temp[:]
            temp=[]
            key_set=[]
    for ky1 in key_set:
        resulting_bins[ky1]=temp[:]
    return resulting_bins


def __bin_withID__(data_dict,min_num=5):
    result=[]
    result_dict={}
    kys=sorted(data_dict.keys())
    temp=[]

    count=0
    for i in range(len(kys)):
	ky=kys[i]
        temp+=data_dict[ky]
        result_dict[ky]=count
        if len(temp)>min_num-1 and i<len(kys)-min_num:
            result.append(temp[:])
            temp=[]
            count+=1
    result.append(temp[:])
    return (result,result_dict)


def comparewithrandom(seedlist,randomsamples,minbin,G,sample_method,sampleArgs=[],statistics=ms.globalmetrics):
    seeds=seedlist[:]

    if sample_method==sampleSnowball:
          temp1=networkx_to_igraph(G,1)
          G_submit=temp1[0]
          no1=temp1[1]
          switch_change=1
          seed_global=statistics(sample_method(G_submit,[no1[x] for x in seeds],sampleArgs))
    else:
          G_submit=G
          switch_change=0
          seed_global=statistics(sample_method(G,seeds,sampleArgs))
    seed_degree=[len(G[x]) for x in seeds] #list of seed's degree
    r=__bin_withID__(__make_num_dict__(G),minbin)
    bins=r[0]
    bins_dict=r[1]
    seed_degree_dict={}
    for item in seed_degree:
        temp_item=bins_dict[item]
        if temp_item in seed_degree_dict:
            seed_degree_dict[temp_item]+=1
        else:
            seed_degree_dict[temp_item]=1
    random_global=[]
    for i in range(randomsamples):
        seeds_temp=sum((rd.sample(bins[y],seed_degree_dict[y]) for y in seed_degree_dict),[])
        if switch_change:
              seeds_temp=[no1[x] for x in seeds_temp]
        random_global.append(statistics(sample_method(G_submit,seeds_temp,sampleArgs)))
    return [seed_global,random_global]



### Compare with Random Wrappers ###

def __multicore_helper__general__(fargs):
	return comparewithrandom(*fargs)

def comparewithrandom_multicore(seedlist,randomsamples,minlist,G,tech,sampleArgs,statistic=ms.globalmetrics,pool1=[]):
    print 'I am assuming i can run this on 20 virtual cores.\nIf this is not the case please cancel this command within 5 seconds.'
    time.sleep(5)
    if pool1==[]:
       p=Pool(processes=20)
    else:
       p=pool1
    f=randomsamples/20
    temp1=[(seedlist,f,minlist,G,sampleArgs) for i in range(20)]
    temp1=[(seedlist,f,minlist,G,functionWrapper(tech),sampleArgs,functionWrapper(statistic)) for i in range(20)]
    r=p.map(__multicore_helper__general__,temp1)
    p.terminate()
    return __glue_and_stick_comparewithrandom__(r)




class functionWrapper():
    def __init__(self,func1):
        self.path='/'.join(x for x in inspect.getmodule(func1).__file__.split('/')[:-1])
        self.module=func1.__module__
        self.name=func1.func_name
        self.f1=None
        self.func_name=self.name
    def __call__(self,*args):
        if not self.f1:
            self.createF1()
        return self.f1(*args)
    def createF1(self):
        sys.path.append(self.path)
        exec('import '+self.module)
        self.f1=eval(self.module+'.'+self.name)



def __glue_and_stick_comparewithrandom__(l1):
	result=l1[0][:]
	for i in range(1,len(l1)):
		result[1]+=l1[i][1]
	return result


###### Networkx to igraph #######

def networkx_to_igraph(G,returnNo=0): #tested
  #  print 'cannot convert to igraph'
    if returnNo==0:
        return G
    else:
        return [G,{x:x for x in G}]
        return (q_G_ig,q_no)
    if G.number_of_nodes()==0:
		if returnNo==0:
			return ig.Graph()
		else:
			return (ig.Graph(),{})

    q_G_ig=ig.Graph()
    q_nodes_d=G.nodes()
    q_no={q_nodes_d[x]:x for x in range(len(q_nodes_d))}
    q_G_ig.add_vertices(max(q_no.values())+1)
    q_G_ig.add_edges((q_no[x[0]],q_no[x[1]]) for x in G.edges())
    q_G_ig.to_undirected()
    if returnNo==0:
        return q_G_ig
    else:
        return (q_G_ig,q_no)



def lcc(G):
    return nx.connected_component_subgraphs(G)[0]
