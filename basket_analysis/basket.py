# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 15:17:40 2015

@author: Wu
"""


import pandas as pd
import numpy as np

import itertools


class NotEnoughData(Exception):
    pass



def kDimNozeroList(NozeroYmasked, k, threshold, ymin = 9, ymax = 21):
    '''
    find the Nozero number in all k-dim combinations
    '''
    result = []
    for index in itertools.combinations(range(ymin, ymax+1), k):
        index_list = list(index)        
        nonzero = np.sum(np.all(NozeroYmasked.loc[:,index_list], 1))
        if nonzero >= threshold:
            result.append((index, nonzero))
    #result.sort(reverse = True, key = lambda x : x[1])
    return result
    #(result, key = lambda x:x[1])
    


def selectAllDim(NozeroY, mask, k, threshold, cluster, cluster_num):
    masked_Y = NozeroY[mask]
    kdim_cluster = kDimNozeroList(masked_Y, k, threshold)
    
    if len(kdim_cluster) == 0:
        raise NotEnoughData
    
    print kdim_cluster
       
    for cluster_i in range(len(kdim_cluster)):
        cluster_num += 1
        
        for data_j in masked_Y.index:
            if all(masked_Y.ix[data_j, list(kdim_cluster[cluster_i][0])]):
                mask[data_j] = False
                cluster[data_j].append(cluster_num)
                
    return kdim_cluster, cluster_num
                



def softCluster(NozeroY, threshold, begin_dim = 7):
    n, m = NozeroY.shape
    
    DataToclusters =[]
    for i in range(n):
        DataToclusters.append([])
    
    clustersMask = np.ones(n, bool)
    cluster_num = 0
    clusters = []
    

    for i in range(begin_dim, 0, -1):
        try:
            kdim_cluster, cluster_num_after = selectAllDim(NozeroY, clustersMask, i ,threshold, DataToclusters, cluster_num)
        except NotEnoughData:
            print("no data enough for dim {0}".format(i))
            continue
        
        for j in range(len(kdim_cluster)):
            index, number = kdim_cluster[j]
            cluster_num += 1
            clusters.append((index, number, cluster_num))
        print('finished dim {0}'.format(i))
           
    return DataToclusters, clusters






if __name__ == "__main__":
    meps = pd.read_csv("2012meps.csv", header = None)
    mepsY = meps.ix[:,9:21]
    NoZeroMepsY = mepsY != 0.0 

    for threshold in [400,700,1000]:
        #threshold = 150
        print('begin {0}'.format(threshold))
        DataTocluster, clusters = softCluster(NoZeroMepsY, threshold)
    
    
        ZeroSum = 0
        NoZeroSum = 0
        for i in range(len(DataTocluster)):
            if len(DataTocluster[i]) == 0:
                if any(NoZeroMepsY.loc[i,:]) != False:
                    NoZeroSum += 1
                    DataTocluster[i].append(-1)
                else:
                    DataTocluster[i].append(0)
                    ZeroSum += 1
    
        
        with open('soft_clustering\\' + str(threshold) +  'DataTocluster.txt', 'w') as fd:
            for i in range(len(DataTocluster)):
                fd.write('{0};{1}\n'.format(i+1, ','.join([str(c) for c in DataTocluster[i]])))
                
        with open('soft_clustering\\' + str(threshold) +  'Clusters.txt', 'w') as fc:
            for c in clusters:
                fc.write('cluster {0}, subspace: {1}, dataInCluster: {2} \n'.format(int(c[2]), c[0], c[1]))
            fc.write('cluster {0}, subspace: {1}, dataInCluster: {2} (allzero)\n'.format(0, 'others', ZeroSum ))
            fc.write('cluster {0}, subspace: {1}, dataInCluster: {2} (notallzero)\n'.format(-1, 'others', NoZeroSum))
        
    
