## LDtools

**LDtools** package for R computes measures and plots of LD for detecting selection. Examples of these tools for LD analysis can be found in the following article:

_Williams F, Meenagh A, Single R, McNally M, Kelly P, Nelson M_P, Meyer D, Lancaster A, Thomson G, and Middleton D. High resolution HLA-DRB1 identification of a Caucasian population. Human Immunology, 2004; 65(1): 66-77._ [http://www.ncbi.nlm.nih.gov/pubmed/14700598](http://www.ncbi.nlm.nih.gov/pubmed/14700598)

**Abstract**
info about DPA ...


### Usage

The **DPAplot** function generates plots and summary measures realted to Disequilibrium Pattern Analysis (DPA), as follows:
```S
dpa.results <- DPAplot(dat=dat, y.threshold=.005, r2.threshold=.70, tolerance=.01)
```

Parameter **dat** is a data.frame with 5 required variables:

  - *haplo.freq* A numeric vector of haplotype frequencies.
  - *locus1* A character vector indentifying the first locus.
  - *locus2* A character vector indentifying the second locus.
  - *allele1* A character vector indentifying the allele at locus 1.
  - *allele2* A character vector indentifying the allele at locus 2.

Parameter **tolerance** is a threshold for the sum of the haplotype frequencies. If the sum of the haplotype frequencies is greater than 1+tolerance or less than 1-tolerance an error is returned.

Parameter **y.threshold** is a threshold for plotting based on the maximum expected freq. If the maximum expected freq is less than y.threshold, no plot is created (default=0.005)

Parameter **r2.threshold** A threshold for plotting based on the fit of the regression line. If the R-squared value is less than r2.threshold, no plot is created (default=0.75)

The function returns a dataframe (in the above example **dpa.results**) with the following components:

- *focal*	 the focal allele at the 1st locus.
- *select*	 the potentially selected allele at the 2nd locus.
- *r2.lt0*	 the R^2 value in the negative D-space.
- *maxdij*	 the maximum d_ij value.
- *exp.frq.max.d*	 the expected freq corresponding to the value with maxdij.
- *prop.gt0*	 the proportion of points with d_ij > 0.
- *n.gt.halfmax.d*  the # of points with d_ij > .5*maxdij.
- *fold.inc*  the fold increase in frequency for the potentially selected haplotype: (hapfreq - expfreq)/expfreq.

### Examples

#### Example 1. HLA frequencies

This is an example of DPA using haplotype frequencies from Williams et al. (2004)

```r
data(NIreland.freqs)
loc1 <- "A"
loc2 <- "B"
temp.dat <- NIreland.freqs[NIreland.freqs$locus1==loc1 & NIreland.freqs$locus2==loc2,]
DPAplot(dat=temp.dat, y.threshold=.005, r2.threshold=.70)
#Create a file with several DPA plots for the chosen loci
postscript(file="Irish_A-B.ps", horizontal=T)
par(mfrow=c(2,2))
DPAplot(dat=temp.dat, y.threshold=.005, r2.threshold=.70)
dev.off()
```

