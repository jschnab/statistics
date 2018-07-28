#script which generates a Mann-Whitney critical U statistics table in which
#the first line is a label containing sample sizes couple [*n1,n2*] followed
#by two lines containing: (1) critical U statistics, (2) their frequencies.
#A third line contains ##### to mark the end of data for this sample sizes
#couple.

import numpy as np
from scipy.stats import rankdata
from scipy.stats import tiecorrect
import itertools

class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1

def MW_gen_tab(x,y):
    ###=== CALCULATION OF U CRITICAL VALUES ===###
    #generates all unique permutations
    fi = open('table_mw2.txt','w')

    #combinations = list(itertools.product(range(x,y),range(x,y)))

    for a,b in combinations:
        #firstly, critical U statistics U_crit and their cumulative frequencies
        #W_cum are calculated
        n_a, n_b = a, b
        liste = list(perm_unique(n_a * [1] + n_b * [0]))

        #transforms unique permutations as list of ranks
        rank_list = np.zeros((len(liste),len(liste[0])))
        j = 0
        while j < len(liste):
            i = 0
            while i < len(liste[j]):
                if liste[j][i]:
                            rank_list[j][i] += i +1
                i +=1
            j +=1

        #calculates W values (sum of ranks) associated with list of ranks
        W_list = np.zeros(len(rank_list))
        j = 0
        while j < len(rank_list):
            W_list[j] = sum(rank_list[j])
            j +=1

        #determines unique W values, their frequency, and cumulated sum of
        #frequencies
        W_unique, W_counts = np.unique(W_list, return_counts =True)
        W_freq = W_counts / len(rank_list)
        W_cum = np.zeros(len(W_freq))
        j = 0
        while j < len(W_freq):
            cum_sum = 0
            for i in range(j +1):
                cum_sum += W_freq[i]
            W_cum[j] = cum_sum
            j +=1

        #calculates critical U statistics
        U_crit = np.zeros(len(W_unique))
        i =0
        while i < len(W_unique):
            U_crit[i] = W_unique[i] - n_a * (n_a +1) / 2
            i +=1

        #secondly, data are written in the file mw_tab.txt
        label = '[*'+str(a)+','+str(b)+'*]'
        fi.write(label +'\n')
        for u in U_crit:
            fi.write(str(u)+',')
        fi.write('\n')
        for f in W_cum:
            fi.write(str(f)+',')
        fi.write('\n')
        fi.write('#####\n')

    fi.close()
