## PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR) model to a collection of training protein sequence alignments.

This analysis fits the GTR model for amino-acids to one or a collection of alignments. The end result is the matrix of (relative  to `I--L`) substitution rates, written in a variety of formats, together with profile likelihood estimates of confidence intervals for said rates.


## Invokation

This analysis has one **required** arguments

- `--list` the list of file paths for the alignments to analyze (paths are relative to the list file, see example below).
 
HyPhy will write Markdown output to the screen and a JSON file with detailed fit results. 


### Complete options list 

Note that `alignment` is only required/used when `--mode Single` is specified. 


```
mode
	Single or Multiple files 
	defaut value: pogofit.multiple [computed at run time]

alignment [required]
	The alignment to analyze

list [required]
	The list of alignments (one per line) to analyze

baseline-model
	The empirical protein model to use for optimizing branch lengths
	defaut value: WAG

frequencies
	Equilibrium frequency estimator
	defaut value: pogofit.emp_freq [computed at run time]

output-format
	Output format for the fitted model
	defaut value: pogofit.output_all [computed at run time]

zero-rates
	Should zero rates be imputed or left at 0
	defaut value: pogofit.impute [computed at run time]

precision
	Optimization precision
	defaut value: 0.001 [computed at run time]
```
 

##Example run 


```
hyphy pogofit_fixgamma.bf --list data/data.list 
```

---
```
Analysis Description
--------------------
PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR)
model to a collection of training protein sequence alignments.

- __Requirements__: All alignments must be in HyPhy-format: Each file must contain a protein
multiple sequence alignment and newick phylogeny. NEXUS input is not
accepted.

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond and Stephanie J Spielman

- __Contact Information__: spond@temple.edu; spielman@rowan.edu

- __Analysis Version__: 0.01

mode: Multiple
list: data/data.list
baseline-model: WAG
frequencies: Emp
output-format: All
zero-rates: Yes


[PHASE 1] Performing initial branch length optimization using WAG
Current Max: -3485.9947     (93 % done) LF Evals/Sec: 2849    CPU Load: 2.638   


[PHASE 2] Optimizing protein model
Current Max: -3379.9298     (99 % done) LF Evals/Sec: 158.4   CPU Load: 1.419   


### Computing confidence intervals for individual rate parameters
A--C:   0.4926 [  0.2121,   0.9553]
A--D:   0.0598 [  0.0000,   0.5739]
A--E:   0.1093 [  0.0084,   0.3082]
A--F:   0.0000 [  0.0000,   0.0846]
A--G:   0.8881 [  0.2226,   2.2668]
A--H:   0.0000 [  0.0000,   0.2342]
A--I:   0.0000 [  0.0000,   0.1677]
A--K:   0.0635 [  0.0000,   0.3103]
A--L:   0.0000 [  0.0000,   0.3181]
A--M:   0.2886 [  0.0615,   0.6990]
A--N:   0.8832 [  0.0000,   4.7966]
A--P:   0.2969 [  0.1214,   0.5739]
A--Q:   0.1816 [  0.0000,   0.8963]
A--R:   0.8233 [  0.1196,   2.4507]
A--S:   0.8449 [  0.4463,   1.4641]
A--T:   0.6448 [  0.0393,   1.7770]
A--V:   0.2189 [  0.0000,   0.9640]
A--W:   0.0993 [  0.0000, 10000.0000]
A--Y:   0.0000 [  0.0000,   7.6353]
C--D:   0.0461 [  0.0000,   0.2487]
C--E:   0.0000 [  0.0000,   0.0470]
C--F:   0.0941 [  0.0122,   0.2512]
C--G:   0.0310 [  0.0000,   0.4472]
C--H:   0.0000 [  0.0000,   0.1866]
C--I:   0.0824 [  0.0000,   0.2423]
C--K:   0.0088 [  0.0000,   0.1240]
C--L:   0.1894 [  0.0000,   0.6311]
C--M:   0.0569 [  0.0000,   0.2706]
C--N:   0.0000 [  0.0000,   1.0831]
C--P:   0.0277 [  0.0000,   0.1059]
C--Q:   0.1158 [  0.0059,   0.4937]
C--R:   0.0000 [  0.0000,   0.4029]
C--S:   0.2938 [  0.1386,   0.5288]
C--T:   0.1076 [  0.0000,   0.5533]
C--V:   0.7202 [  0.2547,   1.5008]
C--W:   0.5906 [  0.0000, 10000.0000]
C--Y:   0.0000 [  0.0000,   7.5149]
D--E:   0.5671 [  0.2959,   1.0120]
D--F:   0.0285 [  0.0000,   0.1826]
D--G:   0.0000 [  0.0000,   0.4598]
D--H:   0.2181 [  0.0000,   0.7764]
D--I:   0.0000 [  0.0000,   0.0767]
D--K:   0.0000 [  0.0000,   0.1733]
D--L:   0.0000 [  0.0000,   0.2797]
D--M:   0.0643 [  0.0025,   0.2932]
D--N:   1.2312 [  0.0363,   6.3015]
D--P:   0.0322 [  0.0000,   0.1937]
D--Q:   0.0000 [  0.0000,   0.6406]
D--R:   1.3426 [  0.3506,   3.4450]
D--S:   0.2419 [  0.0611,   0.5763]
D--T:   0.0025 [  0.0000,   0.5633]
D--V:   0.0000 [  0.0000,   0.3452]
D--W:   0.1137 [  0.0000, 10000.0000]
D--Y:   0.0000 [  0.0000,  12.8304]
E--F:   0.0000 [  0.0000,   0.0343]
E--G:   0.2119 [  0.0398,   0.6145]
E--H:   0.0313 [  0.0000,   0.1788]
E--I:   0.0158 [  0.0006,   0.0687]
E--K:   0.4327 [  0.2615,   0.6846]
E--L:   0.0000 [  0.0000,   0.1226]
E--M:   0.0000 [  0.0000,   0.0480]
E--N:   0.2749 [  0.0000,   1.7967]
E--P:   0.0190 [  0.0007,   0.0670]
E--Q:   0.4515 [  0.1675,   0.9589]
E--R:   0.0000 [  0.0000,   0.2653]
E--S:   0.1048 [  0.0342,   0.2235]
E--T:   0.0000 [  0.0000,   0.1875]
E--V:   0.0000 [  0.0000,   0.1917]
E--W:   0.1368 [  0.0000, 10000.0000]
E--Y:   0.0000 [  0.0000,   2.8677]
F--G:   0.0590 [  0.0000,   0.3463]
F--H:   0.0429 [  0.0000,   0.2017]
F--I:   0.1484 [  0.0386,   0.3289]
F--K:   0.0326 [  0.0000,   0.1280]
F--L:   0.5593 [  0.2020,   1.1922]
F--M:   0.2778 [  0.0988,   0.5783]
F--N:   0.0000 [  0.0000,   0.8410]
F--P:   0.0093 [  0.0000,   0.0516]
F--Q:   0.0000 [  0.0000,   0.1658]
F--R:   0.0000 [  0.0000,   0.2715]
F--S:   0.0353 [  0.0012,   0.1222]
F--T:   0.0000 [  0.0000,   0.1315]
F--V:   0.5348 [  0.1589,   1.2329]
F--W:   1.1628 [  0.0000, 10000.0000]
F--Y:   3.9812 [  0.3730,  45.6237]
G--H:   0.1677 [  0.0005,   0.8203]
G--I:   0.0000 [  0.0000,   0.1679]
G--K:   0.2047 [  0.0264,   0.6858]
G--L:   0.0000 [  0.0000,   0.4971]
G--M:   0.0000 [  0.0000,   0.2403]
G--N:   1.7798 [  0.0352,   9.4779]
G--P:   0.0445 [  0.0000,   0.2776]
G--Q:   0.2650 [  0.0000,   1.4297]
G--R:   0.0000 [  0.0000,   0.9678]
G--S:   0.3837 [  0.0755,   1.0463]
G--T:   0.5016 [  0.0346,   1.7977]
G--V:   0.2039 [  0.0000,   1.3521]
G--W:   0.2888 [  0.0000, 10000.0000]
G--Y:   0.0000 [  0.0000,  24.2199]
H--I:   0.0362 [  0.0000,   0.1726]
H--K:   0.1379 [  0.0000,   0.4592]
H--L:   0.3050 [  0.0157,   0.9429]
H--M:   0.0000 [  0.0000,   0.1238]
H--N:   0.6815 [  0.0000,   4.3691]
H--P:   0.1544 [  0.0402,   0.3619]
H--Q:   0.9821 [  0.3264,   2.2094]
H--R:   1.0861 [  0.2483,   2.9355]
H--S:   0.2846 [  0.1027,   0.6015]
H--T:   0.0000 [  0.0000,   0.2537]
H--V:   0.0000 [  0.0000,   0.3732]
H--W:   0.2269 [  0.0000, 10000.0000]
H--Y:   0.0000 [  0.0000,  10.5900]
I--K:   0.0449 [  0.0000,   0.1610]
I--L: Constrained
I--M:   1.2286 [  0.7793,   1.8780]
I--N:   0.0000 [  0.0000,   0.6673]
I--P:   0.0000 [  0.0000,   0.0304]
I--Q:   0.0000 [  0.0000,   0.1895]
I--R:   0.0000 [  0.0000,   0.2620]
I--S:   0.0325 [  0.0000,   0.1245]
I--T:   0.1015 [  0.0000,   0.4523]
I--V:   2.0123 [  1.1428,   3.3700]
I--W:   0.1846 [  0.0000, 10000.0000]
I--Y:   0.0000 [  0.0000,   5.2901]
K--L:   0.1200 [  0.0000,   0.5013]
K--M:   0.0593 [  0.0000,   0.2201]
K--N:   0.9429 [  0.0960,   3.6228]
K--P:   0.0992 [  0.0349,   0.2067]
K--Q:   0.9238 [  0.3741,   1.8845]
K--R:   1.9035 [  0.9138,   3.6502]
K--S:   0.0409 [  0.0000,   0.1745]
K--T:   0.0000 [  0.0000,   0.2557]
K--V:   0.2953 [  0.0300,   0.8195]
K--W:   0.1204 [  0.0000, 10000.0000]
K--Y:   0.0000 [  0.0000,   5.8013]
L--M:   1.8430 [  0.8831,   3.4313]
L--N:   0.0000 [  0.0000,   2.1383]
L--P:   0.1759 [  0.0522,   0.4115]
L--Q:   0.4006 [  0.0000,   1.6678]
L--R:   0.8185 [  0.0183,   2.6634]
L--S:   0.0000 [  0.0000,   0.1538]
L--T:   0.0000 [  0.0000,   0.5958]
L--V:   1.4027 [  0.3516,   3.5790]
L--W:   0.5509 [  0.0000, 10000.0000]
L--Y:   0.0000 [  0.0000,  15.8193]
M--N:   0.0000 [  0.0000,   1.0910]
M--P:   0.0467 [  0.0014,   0.1430]
M--Q:   0.4673 [  0.1269,   1.1358]
M--R:   0.0000 [  0.0000,   0.5085]
M--S:   0.0000 [  0.0000,   0.1151]
M--T:   0.9838 [  0.4457,   1.8454]
M--V:   0.5202 [  0.0192,   1.6097]
M--W:   0.4337 [  0.0000, 10000.0000]
M--Y:   0.0000 [  0.0000,   8.4175]
N--P:   0.2176 [  0.0000,   1.4225]
N--Q:   0.0000 [  0.0000,   4.3618]
N--R:   0.0000 [  0.0000,   6.8736]
N--S:   0.0000 [  0.0000,   1.1077]
N--T:   0.0000 [  0.0000,   2.8612]
N--V:   0.0000 [  0.0000,   2.8270]
N--W:   0.0635 [  0.0000, 10000.0000]
N--Y:   0.0000 [  0.0000,  64.2488]
P--Q:   0.0000 [  0.0000,   0.1241]
P--R:   0.1943 [  0.0000,   0.6928]
P--S:   0.2480 [  0.1395,   0.4054]
P--T:   0.0000 [  0.0000,   0.2025]
P--V:   0.0200 [  0.0000,   0.1904]
P--W:   0.1331 [  0.0000, 10000.0000]
P--Y:   0.0000 [  0.0000,   3.1076]
Q--R:   1.6749 [  0.2640,   5.0764]
Q--S:   0.0000 [  0.0000,   0.2759]
Q--T:   0.2185 [  0.0000,   1.1901]
Q--V:   0.3338 [  0.0000,   1.5901]
Q--W:   0.1873 [  0.0000, 10000.0000]
Q--Y:   0.0000 [  0.0000,  12.5974]
R--S:   0.1841 [  0.0000,   0.8572]
R--T:   1.0917 [  0.1690,   3.1403]
R--V:   0.0000 [  0.0000,   1.0225]
R--W:   0.9179 [  0.0000, 10000.0000]
R--Y:   0.0000 [  0.0000,  20.3187]
S--T:   1.4639 [  0.8650,   2.3787]
S--V:   0.0000 [  0.0000,   0.2636]
S--W:   0.4392 [  0.0000, 10000.0000]
S--Y:   1.4600 [  0.0568,  14.0603]
T--V:   1.1612 [  0.3219,   2.7802]
T--W:   0.0974 [  0.0000, 10000.0000]
T--Y:   0.0000 [  0.0000,  12.0357]
V--W:   0.3122 [  0.0000, 10000.0000]
V--Y:   0.0000 [  0.0000,  18.6619]
W--Y:   1.7324 [  0.0000, 10000.0000]


 Saving results

```

--- 

##Quick and dirty single file fit run


```
hyphy pogofit_fixgamma.bf --mode Single --alignment data/sim10_treesplit_0.dat --precision 1
```

```
Analysis Description
--------------------
PogoFit, *P*rotein *G*TR *Fit*ter: Fit a general time reversible (GTR)
model to a collection of training protein sequence alignments.

- __Requirements__: All alignments must be in HyPhy-format: Each file must contain a protein
multiple sequence alignment and newick phylogeny. NEXUS input is not
accepted.

- __Citation__: TBD

- __Written by__: Sergei L Kosakovsky Pond and Stephanie J Spielman

- __Contact Information__: spond@temple.edu; spielman@rowan.edu

- __Analysis Version__: 0.01

mode: Single
Provide the filename of the alignment to analyze : alignment: data/sim10_treesplit_0.dat
baseline-model: WAG
frequencies: Emp
output-format: All
zero-rates: Yes
Optimization precision (permissible range = [1e-05,1], default value = 0.001): precision: 1


[PHASE 1] Performing initial branch length optimization using WAG
Current Max: -1401.1448     (98 % done) LF Evals/Sec: 2964    CPU Load: 1.867   


[PHASE 2] Optimizing protein model
Current Max: -1335.2005     (72 % done) LF Evals/Sec: 362.3   CPU Load: 1.506   


### Computing confidence intervals for individual rate parameters
A--C:   1.8137 [  0.2863,   4.5306]
A--D:   0.0000 [  0.0000,   3.5640]
A--E:   0.3840 [  0.0000,   1.6218]
A--F:   0.0000 [  0.0000,   1.0475]
A--G:   2.2219 [  0.0000,   9.4307]
A--H:   0.0000 [  0.0000,   3.7749]
A--I:   0.1975 [  0.0000,   2.3877]
A--K:   0.0000 [  0.0000,   1.8515]
A--L:   0.0000 [  0.0000,   7.2404]
A--M:   1.5869 [  0.0000,   5.8958]
A--N:   0.0000 [  0.0000,  18.0898]
A--P:   1.8483 [  0.4160,   5.0888]
A--Q:   0.0000 [  0.0000,   6.4939]
A--R:   6.0081 [  0.1622,  20.2166]
A--S:   6.5534 [  2.2000,  12.5281]
A--T:   2.9614 [  0.0000,  15.2549]
A--V:   0.0000 [  0.0000,   6.9358]
A--W:   0.0034 [  0.0000, 10000.0000]
A--Y:   0.0000 [  0.0000,  32.0338]
C--D:   0.0000 [  0.0000,   0.9147]
C--E:   0.0000 [  0.0000,   0.3264]
C--F:   0.1405 [  0.0000,   1.0462]
C--G:   0.6639 [  0.0000,   3.4800]
C--H:   0.0735 [  0.0000,   2.0687]
C--I:   0.9777 [  0.1029,   2.3041]
C--K:   0.0000 [  0.0000,   0.6993]
C--L:   0.0000 [  0.0000,   2.8451]
C--M:   0.0000 [  0.0000,   1.4503]
C--N:   0.0000 [  0.0000,   7.4083]
C--P:   0.0000 [  0.0000,   0.3623]
C--Q:   1.7486 [  0.1378,   6.8372]
C--R:   0.0000 [  0.0000,   3.3163]
C--S:   0.6928 [  0.0000,   2.1305]
C--T:   0.3744 [  0.0000,   3.6062]
C--V:   2.6653 [  0.0560,   8.1005]
C--W:   0.0140 [  0.0000, 10000.0000]
C--Y:   0.0000 [  0.0000,  30.3553]
D--E:   0.8901 [  0.1710,   2.5327]
D--F:   0.2770 [  0.0000,   2.0588]
D--G:   0.0000 [  0.0000,   2.9892]
D--H:   2.3368 [  0.0000,  10.4217]
D--I:   0.0000 [  0.0000,   0.8580]
D--K:   0.0000 [  0.0000,   2.3377]
D--L:   0.0000 [  0.0000,   2.7978]
D--M:   0.0000 [  0.0000,   1.4066]
D--N:   0.0000 [  0.0000,  19.1530]
D--P:   0.7623 [  0.0000,   2.9725]
D--Q:   3.2572 [  0.0000,  18.5326]
D--R:   5.2159 [  0.5420,  23.5624]
D--S:   2.0536 [  0.2638,   5.1625]
D--T:   0.8640 [  0.0000,   6.3097]
D--V:   0.0000 [  0.0000,   3.5782]
D--W:   0.0038 [  0.0000, 10000.0000]
D--Y:   0.0000 [  0.0000,  54.9932]
E--F:   0.0000 [  0.0000,   0.2832]
E--G:   0.5344 [  0.0235,   2.4326]
E--H:   0.0000 [  0.0000,   0.8024]
E--I:   0.0000 [  0.0000,   0.2687]
E--K:   1.7269 [  0.9160,   3.3457]
E--L:   0.0000 [  0.0000,   0.8372]
E--M:   0.0000 [  0.0000,   0.4187]
E--N:   1.6697 [  0.0000,  12.8011]
E--P:   0.1385 [  0.0056,   0.5750]
E--Q:   0.0000 [  0.0000,   2.9273]
E--R:   0.0000 [  0.0000,   2.1831]
E--S:   0.4073 [  0.0715,   1.4081]
E--T:   0.1046 [  0.0000,   1.5290]
E--V:   0.2596 [  0.0000,   2.2194]
E--W:   0.0045 [  0.0000, 10000.0000]
E--Y:   0.0000 [  0.0000,  10.0431]
F--G:   0.0000 [  0.0000,   1.5139]
F--H:   0.0000 [  0.0000,   1.7108]
F--I:   0.8622 [  0.1407,   2.7099]
F--K:   0.1956 [  0.0000,   0.9391]
F--L:   0.5734 [  0.0000,   7.2016]
F--M:   1.4867 [  0.3100,   5.6555]
F--N:   0.0000 [  0.0000,   9.8363]
F--P:   0.0000 [  0.0000,   0.4624]
F--Q:   0.0000 [  0.0000,   3.3737]
F--R:   0.0000 [  0.0000,   3.2550]
F--S:   0.0000 [  0.0000,   0.5901]
F--T:   0.0000 [  0.0000,   1.5995]
F--V:   1.7404 [  0.0000,   8.1110]
F--W:   1.7263 [  0.0000, 10000.0000]
F--Y:   6.4543 [  0.2063, 10000.0000]
G--H:   0.0000 [  0.0000,   4.0548]
G--I:   0.0000 [  0.0000,   1.3732]
G--K:   0.0000 [  0.0000,   1.7073]
G--L:   0.0000 [  0.0000,   4.0835]
G--M:   0.0000 [  0.0000,   2.0836]
G--N:   0.0000 [  0.0000,  23.7216]
G--P:   0.4992 [  0.0000,   2.6592]
G--Q:   0.0000 [  0.0000,   7.6006]
G--R:   0.0000 [  0.0000,   7.5426]
G--S:   1.1845 [  0.0000,   5.7100]
G--T:   0.0000 [  0.0000,   4.2615]
G--V:   1.2634 [  0.0000,  10.6705]
G--W:   0.0084 [  0.0000, 10000.0000]
G--Y:   0.0000 [  0.0000,  60.5110]
H--I:   0.0000 [  0.0000,   1.7046]
H--K:   0.1650 [  0.0000,   3.9077]
H--L:   0.0000 [  0.0000,   4.2291]
H--M:   0.0000 [  0.0000,   2.1030]
H--N:   0.0000 [  0.0000,  18.5711]
H--P:   1.4770 [  0.0778,   5.1119]
H--Q:  18.5882 [  5.5919,  55.4540]
H--R:   0.0000 [  0.0000,   9.9957]
H--S:   3.2665 [  1.0995,   8.8739]
H--T:   0.0000 [  0.0000,   4.0992]
H--V:   0.0000 [  0.0000,   5.9056]
H--W:   0.0069 [  0.0000, 10000.0000]
H--Y:   0.0000 [  0.0000,  52.0697]
I--K:   0.3681 [  0.0000,   1.4573]
I--L: Constrained
I--M:  10.0790 [  4.8102,  21.9757]
I--N:   0.0000 [  0.0000,  10.2912]
I--P:   0.0679 [  0.0000,   0.7715]
I--Q:   0.7537 [  0.0000,   5.2119]
I--R:   0.0000 [  0.0000,   5.1310]
I--S:   0.0000 [  0.0000,   0.9873]
I--T:   1.8778 [  0.0000,   6.1775]
I--V:   7.8213 [  2.8599,  28.5623]
I--W:   0.0058 [  0.0000, 10000.0000]
I--Y:   0.0000 [  0.0000,  44.0031]
K--L:   0.0000 [  0.0000,   2.3656]
K--M:   0.0000 [  0.0000,   1.0978]
K--N:   0.0000 [  0.0000,   9.7052]
K--P:   0.4176 [  0.0000,   1.3583]
K--Q:  14.1811 [  5.6629,  39.2034]
K--R:   9.9417 [  3.2653,  27.8947]
K--S:   0.2864 [  0.0000,   1.5830]
K--T:   0.0000 [  0.0000,   1.6007]
K--V:   0.8841 [  0.0000,   5.7976]
K--W:   0.0040 [  0.0000, 10000.0000]
K--Y:   0.0000 [  0.0000,  18.9096]
L--M:  58.9520 [ 16.3465, 270.6814]
L--N:   0.0000 [  0.0000,  39.5617]
L--P:   0.9972 [  0.0779,   4.9460]
L--Q:   0.0000 [  0.0000,  10.2162]
L--R:   8.7277 [  1.5526,  42.7262]
L--S:   0.0000 [  0.0000,   2.0734]
L--T:   0.0000 [  0.0000,  13.7570]
L--V:   5.0226 [  0.0000,  53.6161]
L--W:   0.0134 [  0.0000, 10000.0000]
L--Y:   0.0000 [  0.0000, 241.4394]
M--N:   5.9593 [  0.0000,  29.7052]
M--P:   0.2920 [  0.0000,   1.8546]
M--Q:   0.0000 [  0.0000,   4.7507]
M--R:   0.8058 [  0.0000,  13.8808]
M--S:   0.0000 [  0.0000,   1.1242]
M--T:   3.6153 [  0.0000,  11.9295]
M--V:   4.6668 [  0.0000,  25.9198]
M--W:   0.0141 [  0.0000, 10000.0000]
M--Y:   0.0000 [  0.0000, 137.6139]
N--P:   0.0000 [  0.0000,   6.9019]
N--Q:   0.0000 [  0.0000,  43.1101]
N--R:   0.0000 [  0.0000,  60.7334]
N--S:   0.0000 [  0.0000,   9.5075]
N--T:   6.7533 [  0.0000,  62.5554]
N--V:   0.0000 [  0.0000,  38.5628]
N--W:   0.0022 [  0.0000, 10000.0000]
N--Y:   0.0000 [  0.0000, 478.6291]
P--Q:   0.0000 [  0.0000,   3.8783]
P--R:   0.8019 [  0.0000,   6.6505]
P--S:   0.7610 [  0.1305,   2.6100]
P--T:   0.9505 [  0.0000,   3.8449]
P--V:   0.0000 [  0.0000,   1.7965]
P--W:   0.0044 [  0.0000, 10000.0000]
P--Y:   0.0000 [  0.0000,  14.8087]
Q--R:   7.6867 [  0.0000,  53.7971]
Q--S:   0.0000 [  0.0000,   5.6289]
Q--T:   0.0000 [  0.0000,   7.1226]
Q--V:   0.0000 [  0.0000,  14.9436]
Q--W:   0.0059 [  0.0000, 10000.0000]
Q--Y:   0.0000 [  0.0000,  90.8210]
R--S:   0.0000 [  0.0000,   4.3314]
R--T:   5.9849 [  0.0000,  22.2372]
R--V:   0.0000 [  0.0000,  17.5584]
R--W:   1.7371 [  0.0000, 10000.0000]
R--Y:   0.0000 [  0.0000,  86.0668]
S--T:   5.1383 [  1.9614,  13.6159]
S--V:   0.8939 [  0.0000,   5.8143]
S--W:   0.0115 [  0.0000, 10000.0000]
S--Y:   0.0000 [  0.0000,  17.7475]
T--V:   5.4299 [  0.4952,  29.8848]
T--W:   0.0033 [  0.0000, 10000.0000]
T--Y:   0.0000 [  0.0000,  46.9415]
V--W:   0.0089 [  0.0000, 10000.0000]
V--Y:   0.0000 [  0.0000, 166.6773]
W--Y:   1.7089 [  0.0000, 10000.0000]


 Saving results


Analysis complete!
```
