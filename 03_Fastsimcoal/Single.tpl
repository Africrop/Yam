//Number of population samples (demes)
4

//Population effective sizes (number of genes)
NPOP0$
NPOP1$
NPOP2$
NPOP3$ 

//Sample sizes
20
20
20
20

//Growth rates: negative growth implies population expansion
0
0
0
0

//Number of migration matrices : 0 implies no migration between demes
0

//Historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
// Model 01
3  historical event 
T0$ 3 0 1 1 0 0
T1$ 1 2 1 RES1$ 0 0   
T2$ 2 0 1 RES2$ 0 0  

// Model 02
//3  historical event 
//T0$ 3 1 1 1 0 0
//T1$ 1 2 1 RES1$ 0 0   
//T2$ 2 0 1 RES2$ 0 0 

//Model 03
//3  historical event 
//T0$ 3 2 1 1 0 0
//T1$ 1 2 1 RES1$ 0 0   
//T2$ 2 0 1 RES2$ 0 0 

//Number of independent loci [chromosome]
1 0

//Per chromosome: Number of linkage blocks
1

//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 MutRate$
