# mstree

This Perl program calculate the time to the most recent common ancestor (*TMRCA*) and minimum clade proportion (*p<sub>mc</sub>*) for coalescent genealogies generated with [ms](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html).

These values were reported in:

Cox MP, Mendez FL, Karafet TM, Metni Pilkington M, Kingan SB, Destro-Bisol G, Strassmann BI, Hammer MF. 2008. [Testing for archaic hominin admixture on the X chromosome: Model likelihoods for the modern human *RRM2P4* region from summaries of genealogical topology under the structured coalescent](https://doi.org/10.1534/genetics.107.080432). *Genetics* 178: 427-437.

Specifically, the code implements equation A1 in Cox *et al* (2008) for the minimum clade proportion (*p<sub>mc</sub>*):

<img src="Cox_EquationA1.jpg" width="300"/>







The [*gst_prime* function](gst_prime.R) is written in base R, and requires values for the number of subpopulations *k*, the average subpopulation heterozygosity *H<sub>S</sub>* and the total population heterozygosity *H<sub>T</sub>*.  Usage is:

gst_prime(*k*, *H<sub>S</sub>*, *H<sub>T</sub>*)

Worked example:

A study population containing 14 subpopulations with average subpopulation heterozygosity of 0.953 and total population heterozygosity of 0.981 would return:

```
gst_prime(14, 0.953, 0.981)
0.6518016
```