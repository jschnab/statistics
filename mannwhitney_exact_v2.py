
# 01/03/2018
# this module calculates the exact p-value from a two-sided Mann-Whitney test
# (no longer the case)the Mann-Whitney test from scipy.stats is used to
# calculate the U statistics
# the U statistics is then compared to a Wilcoxon distribution
# (no longer the case) the distribution is calculated for each test, so samples
# should be kept small, i.e. < 10

# 04/03/2018
# adding code which reads Mann-Whitney U and p-values from a text file called
# table_mw.txt

# 07/03/2018
# removed un-necessary instructions previously used to generate U statistics
# tables.

# 29/07/2018
# removed note about critical values calculated on the fly (they are now read from
# a file.

import numpy as np
from scipy.stats import rankdata
from scipy.stats import tiecorrect

def MW_exact(a, b, alternative ='two-sided'):
    """
    Performs a two-sided Mann-Whitney test with two samples.

    Parameters
    ----------
    a, b : list
        List of ordinal values corresponding to two samples.
    alternative : 'less', 'two-sided', or 'greater'
        Whether to calculate the p-value for the one-sided hypothesis
        ('less' or 'greater') or for the two-sided hypothesis ('two-sided').
        Default to 'two-sided'.

    Returns
    -------
    U statistic : float
        U statistic from the Mann-Whitney test, equal to U for sample a.
    p-value : float
        Returns the exact p-value. If there are ties the indicated p-value corresponds
        to the U statistic rounded to the next integer (i.e. 5 if U = 4.5).
    U_crit : array of float
        Critical U values calculated given the sample sizes.
    W_cum : array of float
        Cumulative frequency of Wilcoxon W statistic (or U statistic) given the
        sample sizes.

    Notes
    -----

    This function calculates only exact values, so if there are ties the correct
    U statistics is returned but p-value corresponds to U rounded to the next integer
    (i.e. 5 if U = 4.5).

    References
    ----------
    .. [1] http://www.real-statistics.com/non-parametric-tests/wilcoxon-rank-sum-test/wilcoxon-rank-sum-exact-test/
    
    .. [2] http://www.real-statistics.com/non-parametric-tests/mann-whitney-test/mann-whitney-exact-test/

    .. [3] https://en.wikipedia.org/wiki/Mann-Whitney_U_test

    .. [4] H.B. Mann and D.R. Whitney, "On a Test of Whether one of Two Random
           Variables is Stochastically Larger than the Other," The Annals of
           Mathematical Statistics, vol. 18, no. 1, pp. 50-60, 1947.
    """
    ###=== READING OF U statistics and p-values FROM TEXT FILE ===###
    #opens a text file containing Mann-Whitney tables and makes
    #a dictionary dic_tot in which keys are sample size couples and items are
    #dictionaries in which keys are U_crit and W_cum and items the corresponding
    #values
    fi=open('table_mw.txt','r')

    dic_tot = {} #dico_tot['3,3']['U_crit']
    dic = {}     #dic['U_crit'] and dic['W_cum']
    lines = []
    i = 0

    for line in fi:
        if line[:2]=='[*':
            label=line[2:]
            label=label[:-1].strip()
            label=label[:-2]
        else:
            if line[:5]=='#####':
                dic['U_crit'] = lines[0]
                dic['W_cum'] = lines[1]
                dic_tot[label]=dic
                lines = []
                i = 0
                dic = {}
            else:
                lines.append(line)
    fi.close()

    #extracts U_crit and W_cum values from dictionaries
    n_a, n_b = len(a), len(b)
    U_crit_list = dic_tot[str(n_a)+','+str(n_b)]['U_crit'].split(',')
    U_crit = [float(i) for i in U_crit_list[:-1]]
    W_cum_list = dic_tot[str(n_a)+','+str(n_b)]['W_cum'].split(',')
    W_cum = [float(i) for i in W_cum_list[:-1]]

    ###=== MANN-WHITNEY TEST ===###
    #calculates U statistics from Mann-Whitney test then compares it to
    #critical values
    a = np.asarray(a)
    b = np.asarray(b)
    ranked = rankdata(np.concatenate((a,b)))
    rank_a = ranked[:n_a] #get ranks for sample a
    u_a = np.sum(rank_a, axis =0) -n_a*(n_a +1)/2
    u_b = n_a*n_b -u_a #remember u_a +u_b == n_a*n_b

    if alternative == 'two-sided':
        U = min(u_a, u_b)
    elif alternative == 'less':
        U = u_a
    elif alternative == 'greater':
        U = u_b

    if tiecorrect(ranked) == 0:
        raise ValueError('All numbers are identical')
    
    i = 0
    while i < len(U_crit):
        p_value = W_cum[i]
        if U <= U_crit[i]:
            break
        i +=1

    if alternative == 'two-sided':
        p_value = p_value*2
    
    return U, p_value, U_crit, W_cum

if __name__ == '__main__':
    from random import randint
    a ,b = [], []
    for i in range(7):
        a.append(randint(20,80))
    for i in range(7):
        b.append(randint(40,100))
    print('Samples : ',a,b)
    print('Mann-Whitney two-sided : ',end='')
    print(MW_exact(a,b)[:2])
    print('Mann-Whitney one-sided (less) : ',end='')
    print(MW_exact(a,b,alternative='less')[:2])
    print('Mann-Whitney one-sided (greater) : ',end='')
    print(MW_exact(a,b,alternative='greater')[:2])
