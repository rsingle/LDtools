## LDtools

**LDtools** package for R computes measures and plots of LD for detecting selection. Examples of these tools for LD analysis can be found in the following article:

_Williams F, Meenagh A, Single R, McNally M, Kelly P, Nelson M_P, Meyer D, Lancaster A, Thomson G, and Middleton D. High resolution HLA-DRB1 identification of a Caucasian population. Human Immunology, 2004; 65(1): 66-77._ [http://www.ncbi.nlm.nih.gov/pubmed/14700598](http://www.ncbi.nlm.nih.gov/pubmed/14700598)

**Abstract**
info about DPA ...


### Usage

To compute ALD in a data set, you can use the function **compute.ALD** within the asymLD package, as follows:
```S
ald.results <- compute.ALD(dat, tolerance = 0.005)
```

Parameter **dat** is a data.frame with 5 required variables:

  - *haplo.freq* A numeric vector of haplotype frequencies.
  - *locus1* A character vector indentifying the first locus.
  - *locus2* A character vector indentifying the second locus.
  - *allele1* A character vector indentifying the allele at locus 1.
  - *allele2* A character vector indentifying the allele at locus 2.

Parameter **tolerance** is a threshold for the sum of the haplotype frequencies. If the sum of the haplotype frequencies is greater than 1+tolerance or less than 1-tolerance an error is returned.

Parameter **symm** is an indicator for whether to compute symmetric measures Dprime & Wn (default=FALSE)

The function returns a dataframe (in the above example **ald.results**) with the following components:

- *locus1*	The name of the first locus.
- *locus2*	The name of the second locus.
- *F.1*	Homozygosity (expected under HWP) for locus 1.
- *F.1.2*	Conditional homozygosity for locus1 given locus2.
- *F.2*	Homozygosity (expected under HWP) for locus 2.
- *F.2.1*	Conditional homozygosity for locus2 given locus1.
- *ALD.1.2*	Asymmetric LD for locus1 given locus2.
- *ALD.2.1*	Asymmetric LD for locus2 given locus1.


### Examples

#### Example 1. HLA frequencies

This is an example of ALD measure using haplotype frequencies from Wilson (2010)

```r
library(asymLD)
data(hla.freqs)
hla.a_b <- hla.freqs[hla.freqs$locus1=="A" & hla.freqs$locus2=="B",]
head(hla.a_b)
```

```
##     haplo.freq locus1 locus2 allele1 allele2
## 170    0.06189      A      B    0101    0801
## 171    0.04563      A      B    0201    4402
## 172    0.04318      A      B    0301    0702
## 173    0.03103      A      B    0201    4001
## 174    0.02761      A      B    0301    3501
## 175    0.01929      A      B    0205    1503
```

```r
compute.ALD(hla.a_b)
```

```
##   locus1 locus2       F.1    F.1.2        F.2     F.2.1   ALD.1.2  ALD.2.1
## 1      A      B 0.1021811 0.340332 0.04876543 0.1903394 0.5150291 0.3857872
```

```r
hla.c_b <- hla.freqs[hla.freqs$locus1=="C" & hla.freqs$locus2=="B",]
compute.ALD(hla.c_b)
```
```
##   locus1 locus2        F.1     F.1.2        F.2     F.2.1   ALD.1.2    ALD.2.1
## 1      C      B 0.08637241 0.7350216 0.05040254 0.4520268 0.8425979  0.6503396
```

Note that there is substantially less variablity (higher ALD) for HLA\*C conditional on HLA\*B than for HLA\*B conditional on HLA\*C, indicating that the overall variation for C is relatively low given specific B alleles.
