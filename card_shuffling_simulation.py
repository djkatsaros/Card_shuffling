#! /usr/bin/python -w

# Card Shuffling "cut-off" phenomenon demonstration
# 
# Based on the discussion in:
#   https://statweb.stanford.edu/~cgates/PERSI/papers/aldous86.pdf
#   of the question "how many shuffles does it take to be 'sufficiently' random?"
#
# 'Sufficiently random' is tken to mean that we are close to a uniform 
#   distribution on the cards in the deck.
#
# The script demonstrates that a sharp bound exists for a given deck size. I.e., 
#   for 10 cards, 5 shuffles might not get us very close to uniform, but upon shuffling
#   a 6th time, we are very close to uniform. This is the so called "cut-off"
#   phenomenon. 
#
# The following script calculates and plots the total variational difference 
#   from uniform aftter  k 'riffle' shuffles [traditional shuffling technique where one
#   splits the deck into two halves and then 'riffles' them together].
#
# It turns out that the probability that we move from the indentity to a given permutation
#   is a function of the number of rising sequences of that permutation.
#   see: https://pdfs.semanticscholar.org/24dc/3948bd857a113ba4fa80cd241f6f60c3133a.pdf
#   Thus, the script incldues a method that computes the number of rising sequences of a 
#   given permutation.
#   While this method is computationally expensive, we note that this is less cumbersome th#   having to construct  a proper random walk on a k! state space - k the size of the deck.
#
#   For a further exposition of this topic and the project that led to this script, see
#       https://people.umass.edu/dkatsaros/card_shuffling_analysis.pdf


import numpy as np
import math
import itertools as it
from matplotlib import pyplot as plt


def binom(n,k):
    """ Computes the binomial coefficient for intergers $n, k$. If you don't know
    the binom. coeff., idk why you're here, but just in case:
    https://en.wikipedia.org/wiki/Binomial_coefficient
        """
    if n == k:
        out = 1
    elif k == 1:
        out = n
    elif k > n:
        out = 0
    else:
        a = math.factorial(n)
        b = math.factorial(k)
        c = math.factorial(n-k)
        div = a // (b * c)
        out = div
    return out

def prob_shuff(sig, k, n):
    """ Probabilty of going from the null permutation to a give permutation in k
    2-shuffles (Corresponding to splitting deck into 2 halves and riffle shuffling, GSR model) 
    > n the number of cards in a deck labeled 1 thru n.
    > sig is a permutation stored as a list
    > k number of shuffles
   """
    r = rising_seqs(sig, {1 : [sig[0]]})
    Qk = 2**(-k*n) * binom(n +2**k - r, n)
    return Qk


def rising_seqs(l,rs):
    """ Computes and enumerates the number of rising sequences of a list  l
        See: https://statweb.stanford.edu/~cgates/PERSI/papers/bayer92.pdf
        for definition and examples of rising sequencesi
        > l a permutation, inputed as a list
        > rs a dictionary 
        """
    for i in np.arange(1,np.shape(l)[0]):
        cand = l[i]
        ### Uncomment to enumerate the candidate sequences -> will dominate your console
        #print('candidate sequence=',cand)
        new = True
        for k in list(rs.keys()):
            #print('k=',k)
            s = rs[k]
            if cand - s[-1] == 1:
                s = list(np.append(s,cand))
                #print('s=',s)
                rs[k] = s
                new = False
        if new == True:
            rs[max(list(rs.keys())) +1] = [cand]

        out = len(rs)
    return out

def total_variation(l1,l2):
    """ takes in two "distributions" (lists l1 l2) and returns the total total_variation
        distance  between them"""
    #if np.size(l1) == np.size(l2)
    nfact = np.size(l1)
    diffs = [0] * nfact
    for j in np.arange(0,nfact):
        diffs[j] = abs(l1[j] - l2[j])

    TV = 0.5* sum(diffs)

    return TV


def diff(n,k):
    """ calculates the total variation difference from uniform after k shuffles.
        Not vectorized. 
    > input size of deck and number of shuffles
    > n = int(input("Enter number of cards: "))
    > k = int(input("Enter how many times shuffle: "))
    > null permutation """

    # initial permuation: just 1, 2, \ldots, n 
    sig0 = list(np.arange(1,n+1))
    # generates all permuations
    Sn = list(it.permutations(sig0))
    # stores numbers
    nfact  = math.factorial(n)
    U = [1/nfact] * nfact
    distr = [0] * nfact
 
    for i  in np.arange(0, nfact):
        sig = Sn[i]
        distr[i] = prob_shuff(sig, int(k), n)

    diff = total_variation(distr, U)

    return diff,distr

### Some scripted inputs, if you're into that
#err = 3*1e-3
#n = 4

def main(n,err):
    k = 1
    diffr, distr = diff(n,k)
    distrs = [distr[0:9]]
    diffrs = [diffr]

    while diffr >= err:
        print('for ', k, ' shuffles, get total variation of %1.7f , [out of tolerance] ' % diffr)
        diffr, distr= diff(n,k)
        diffrs = np.append(diffrs,diffr)
        distrs = np.append(distrs, [distr[0:9]])
        k = k+1

    print('for ', k, ' shuffles, get total variation of %1.7f , Stop, within tolerance.' % (diffr))
    distrs = list(np.around(np.array(distrs),8))
    for k in [0,1,2,3,4,5]:
        print(list(distrs[0+10*k:9+10*k]))
    
    # plotting, to see cut-off pheonomenon
    plt.plot(diffrs)
    plt.show()
    return distrs


if __name__ == "__main__":
    distrs = main(n,err)



#-------------------------------------------------------------------------------
