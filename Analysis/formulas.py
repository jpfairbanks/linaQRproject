from __future__ import print_function
import itertools
from itertools import starmap
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as plt
import sys
def spgivenscostvec(m,n,i):
	"""Compute the cost of using givens rotations
	to eliminate the ith column of m,n matrix"""
	elimcost = 6.0*(m-i)*(n-i)
	qcost = 12.0*(m-i)*n
	fillcost = 3.0*(n-i)*((n-i)-1)
	matveccost = m*m
	return elimcost,qcost,fillcost, #matveccost

def allpairssimulate(ms, ns, iset):
	for m in ms:
		for n in ns:
			for i in isset:
				if i<=n:
					hc = sum(householdercost(m,n))
					gc = sum(spgivenscostvec(m,n,1))
					yield((m,n,i, i*gc/hc))

def geninput():
	ms = [100*x for x in range(10)]
	ns = [20*x for x in range(10)]
	i = [n/2 for n in ns]
	print(ms,ns,i)
	return ms,ns,i

def householdercost(m,n):
	accumcost = 4*(m*m*n - m*n*n + n*n*n/3)
	maincost  = (2*(n*n)) + (m-n)**3
	return maincost, accumcost

def printl(seq):
	for item in seq:
		print (item)
	return seq
params = [(100*x, 20*x, 10*x) for x in np.arange(10)+1]
print("m,n,i")
printl(params)
costs = [spgivenscostvec(p[0], p[1], p[2]) for p in params]
hhcosts = [householdercost(pt[0], pt[1]) for pt in params]
diffs = [(sum(hc)/sum(gc)) for gc,hc in zip(costs, hhcosts)]
print("Givens")
[print((c, sum(c))) for c in costs]
print("Householder")
[print(c) for c in hhcosts]
print("HH-GR/GR")
[print(d) for d in diffs]

ms = [100, 500, 1000]
ns = [16* 2**i for i in range(6)]
isset = [2**i for i in range(10)]

printl(allpairssimulate(ms, ns, isset))
