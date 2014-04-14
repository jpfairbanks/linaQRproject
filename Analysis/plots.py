from __future__ import print_function
import itertools
from itertools import starmap
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import sys

def speedupframe(df, m):
	"""For a fixed value of m, create a DataFrame indexed by
	values of k where each column represents a value of n
	and the entries are the speedup ratios."""
	frame = df[df.m == m]
	sframe = frame[["m","n","k","speedup"]]
	sframe.index = sframe.k
	#print(sframe)
	#sframe[sframe.n == n][["k","speedup"]]
	#"n={0}".format(n)
	nseqs = {n:sframe[sframe.n == n].speedup
	 			for n in np.unique(sframe.n)}
	ratioframe = pd.DataFrame(nseqs)
	#print(ratioframe)
	return ratioframe


def plotspeedup(spframe, m):
	"""Plot the frame with appropriate titles and axes"""
	spframe["breakeven"] = 1
	ax = spframe.plot(logy=False, logx=False)
	ax.set_ylabel("Time Full / Time Eager")
	ax.set_title("Runtime Ratios m = {0}".format(m))
	return ax

def modeladjustment(sf, m):
	""" Adjust the speedup factors to account for the
	analytic flop counts upto bigO. 

	tfull/teager = C * n/k

	C = tfull*k  /  teager*n 

	Returns a DataFrame indexed by k with n values for columns.
	Each entry is the estimate of C by solving the above equation
	at the value of n,k. 
	"""
	dd = {n:sf[n] * sf[n].index / (int(n)+0.0) for n in sf.columns}
	adjf = pd.DataFrame(dd)
	return adjf


def predictBreakeven(sf, m):
	"""Use the analytic flops counts bigO to predict
	 the break even point based on speedups.
	 sf should be retured by speedupframe

	 Returns a Series indexed by n with entries of k
	"""
	adjf = modeladjustment(sf, m)
	Cseq = adjf.mean()
	predbreakeven = Cseq*Cseq.index
	return predbreakeven

def breakevenpt(df,m,n):
	""" Provide bounds on k_even := smallest k where tratio(m,n,k) < 1
	example: lb, ub = breakevenpt(df, 500, 32)
	then lb <= k_even <= ub

	This is entirely empirical. You can check that the
	estimates of breakevenpt are within the bounds by 
	using predictBreakeven to get model aware estimates of
	the break even points.
	"""
	mnframe = df[df.m == m][df.n == n][["k","speedup"]]
	lossframe = mnframe[mnframe.speedup <= 1.0]
	gainframe = mnframe[mnframe.speedup >= 1.0]
	ub = lossframe.k.min()
	lb = gainframe.k.max()
	return lb, ub
def computeBounds(df):
	"""Applies breakevenpt to find an upper bound
	and lower bound for the break even point value of k
	for all pairs of m,n values.

	Because of the choices of square cases and 
	discretization from powers of 2 not all cases will be run.
	"""
	mlevels = np.unique(df.m)
	nlevels = np.unique(df.n)
	#allocate
	lb = np.zeros((len(mlevels), len(nlevels)))
	ub = np.zeros((len(mlevels), len(nlevels)))
	for i,m in enumerate(mlevels):
		for j,n in enumerate(nlevels):
			lb[i,j], ub[i,j] = breakevenpt(df,m,n)
	return lb, ub

if __name__ == '__main__':
	#fp = "datafull.tsv"
	fp = "workhorseoutput.tsv"
	plots = True
	breakeven = True
	model = True
	#colnames = ["m", "n", "k", "tfull","teager","tlazy","sigfull","sigeager","siglazy","pfe", "pel", "pfl"]
	df = pd.read_csv(fp, 
		header=0,# names=colnames,
		index_col=False,
		sep="\t",)
	#df = pd.read_csv(fp, header=0, index_col=0)
	mlevels = np.unique(df.m)
	nlevels = np.unique(df.n)
	print(df)
	df["speedup"] = df.tfull/df.teager
	if plots:
		for m in mlevels:
			sf = speedupframe(df,m)
			ax = plotspeedup(sf,m)
			ax.figure.savefig("tratio{0}.png".format(m))
	if model:
		kevendict = {}
		for m in mlevels:
			sf = speedupframe(df,m)
			kevendict[m] = predictBreakeven(sf,m)
			# print("The break even points for m={0}".format(m))
			# print(keven)
		kevenf = pd.DataFrame(kevendict)
		print(kevenf)
	if breakeven:
		lb,ub = computeBounds(df)
		print(lb.T)
		print(ub.T)
		mids = (lb+ub)/2
		mids = mids.T
		print(mids)