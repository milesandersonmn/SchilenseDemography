//Number of population samples (demes)
7 samples to simulate :
//Population effective sizes (number of genes)
NCV1$
NCV2$
NCV3$
NCV4$
NSCG$
NSHG$
NGHOST$
//Samples sizes
20
18
18
20
12
18
0
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
6 historical event
TDIV1$ 0 6 1 RES1$ 0 0
TDIV2$ 1 6 1 RES2$ 0 0
TDIV3$ 2 6 1 RES3$ 0 0
TDIV4$ 3 6 1 RES4$ 0 0
TDIVSHG$ 5 6 1 RES5$ 0 0
TDIVSCG$ 4 6 1 RES6$ 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0.0000 0.00000001 OUTEXP
