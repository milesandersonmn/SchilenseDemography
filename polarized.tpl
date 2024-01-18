//Number of population samples (demes)
38 samples to simulate :
//Population effective sizes (number of genes)
NLA1958$ sfspool 0
NLA1960$ sfspool 0
NLA2981A$ sfspool 0
NLA3115$ sfspool 1
NLA1968$ sfspool 1
NLA3114$ sfspool 1
NLA3112$ sfspool 1
NLA3111$ sfspool 1
NLA1963$ sfspool 2
NLA1967$ sfspool 2
NLA2965$ sfspool 2
NLA1965$ sfspool 2
NLA2751$ sfspool 3
NLA2753$ sfspool 3
NLA2754$ sfspool 3
NLA2755$ sfspool 3
NLA4109$ sfspool 4
NLA4108$ sfspool 4
NLA4107$ sfspool 4
NLA4106$ sfspool 4
NLA4338$ sfspool 4
NLA4339$ sfspool 4
NLA2884$ sfspool 5
NLA4329$ sfspool 5
NLA4330$ sfspool 5
NLA4332$ sfspool 5
NLA4117A$ sfspool 5
NLA4118$ sfspool 5
NLA4119$ sfspool 5
NLA4335$ sfspool 5
NLA2880$ sfspool 5
NCV1$ sfspool -1
NCV2$ sfspool -1
NCV3$ sfspool -1
NCV4$ sfspool -1
NSCG$ sfspool -1
NSHG$ sfspool -1
NGHOST$ sfspool -1
//Samples sizes
20
18
18
18
20
18
20
20
18
16
18
18
20
18
18
20
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
0
0
0
0
0
0
0
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
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
TJOINCV1$ 0 31 1 1 0 0
TJOINCV1$ 1 31 1 1 0 0
TJOINCV1$ 2 31 1 1 0 0
TJOINCV2$ 3 32 1 1 0 0
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
