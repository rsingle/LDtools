#' Create DPA plots.
#'
#' A function to generate Disequilibrium Patern Analysis (DPA) plots for haplotype frequency data. 
#' 
#' @param dat	A data.frame with 5 required variables (having the names listed below):
#'      \tabular{ll}{
#'           \code{haplo.freq} \tab A numeric vector of haplotype frequencies.\cr
#'            \code{locus1} \tab A character vector indentifying the first locus.\cr
#'            \code{locus2} \tab A character vector indentifying the second locus.\cr
#'            \code{allele1} \tab A character vector indentifying the allele at locus 1.\cr
#'            \code{allele2} \tab A character vector indentifying the allele at locus 2.\cr
#'            }
#' @param tolerance A threshold for the sum of the haplotype frequencies.
#'    If the sum of the haplotype frequencies is greater than 1+tolerance or less
#'    than 1-tolerance an error is returned. The default is 0.01.
#'
#' @param y.threshold A threshold for plotting based on the maximum expected freq. 
#'    If the maximum expected freq is less than y.threshold, no plot is created (default=0.005)
#'
#' @param r2.threshold A threshold for plotting based on the fit of the regression line.
#'    If the R-squared value is less than r2.threshold, no plot is created (default=0.75)

#' @return A series of plots are created. The return value is a dataframe with the following components:
#'  \tabular{ll}{
#'  \code{focal}	\tab the focal allele at the 1st locus.\cr
#'  \code{select}	\tab the potentially selected allele at the 2nd locus.\cr
#'  \code{r2.lt0}	\tab the R^2 value in the negative D-space.\cr
#'  \code{maxdij}	\tab the maximum d_ij value.\cr
#'  \code{exp.frq.max.d}	\tab the expected freq corresponding to the value with maxdij.\cr
#'  \code{prop.gt0}	\tab the proportion of points with d_ij > 0.\cr
#'  \code{n.gt.halfmax.d} \tab the # of points with d_ij > .5*maxdij.\cr
#'  \code{fold.inc} \tab the fold increase in frequency for the potentially selected haplotype: (hapfreq - expfreq)/expfreq.\cr
#'  }
#'
#' @section Details:
#' A warning message is given if the sum of the haplotype frequencies is greater than 1.01 or less
#' than 0.99 (regardless of the \code{tolerance} setting). The haplotype frequencies that are passed
#' to the function are normalized within the function to sum to 1.0 by dividing each frequency by 
#' the sum of the passed frequencies.
#'
#' @examples
#' library(LDtools)
#' 
#' # An example using the Northern Ireland data from Williams et al.(2004)
#' data(NIreland.freqs)
#' loc1 <- "A"
#' loc2 <- "B"
#' temp.dat <- ni.dat[ni.dat$locus1==loc1 & ni.dat$locus2==loc2,]
#' DPAplot(dat=temp.dat, y.threshold=.005, r2.threshold=.70)
#' #Create a file with several DPA plots for the chosen loci
#' postscript(file="Irish_A-B.ps", horizontal=T)
#' par(mfrow=c(2,2))
#' DPAplot(dat=temp.dat, y.threshold=.005, r2.threshold=.70)
#' dev.off()
#'
#' #' # An example using haplotype frequencies from Wilson(2010)
#' data(hla.freqs)
#' hla.a_b <- hla.freqs[hla.freqs$locus1=="A" & hla.freqs$locus2=="B",]
#' compute.ALD(hla.a_b)
#' hla.freqs$locus <- paste(hla.freqs$locus1, hla.freqs$locus2, sep="-")
#' compute.ALD(hla.freqs[hla.freqs$locus=="C-B",])
#' # Note: additonal columns on the input dataframe (e.g., "locus" above) are allowed, but 
#' # ignored by the function.
#' 
#' # An example using genotype data from the haplo.stats package
#' require(haplo.stats)
#' data(hla.demo)
#' geno <- hla.demo[,5:8]  #DPB-DPA 
#' label <- unique(gsub(".a(1|2)", "", colnames(geno)))
#' label <- paste("HLA*",label,sep="")
#' keep <- !apply(is.na(geno) | geno==0, 1, any)
#' em.keep  <- haplo.em(geno=geno[keep,], locus.label=label)
#' hapfreqs.df <- cbind(em.keep$haplotype, em.keep$hap.prob) 
#' #format dataframe for ALD function
#' names(hapfreqs.df)[dim(hapfreqs.df)[2]] <- "haplo.freq"
#' names(hapfreqs.df)[1] <- "allele1"
#' names(hapfreqs.df)[2] <- "allele2"
#' hapfreqs.df$locus1 <- label[1]
#' hapfreqs.df$locus2 <- label[2]
#' head(hapfreqs.df)
#' compute.ALD(hapfreqs.df)
#' # Note that there is substantially less variablity (higher ALD) for HLA*DPA1 
#' # conditional on HLA*DPB1 than for HLA*DPB1 conditional on HLA*DPA1, indicating 
#' # that the overall variation for DPA1 is relatively low given specific DPB1 alleles
#' 
#' @export
#' @importFrom stats aggregate

DPAplot <- function(dat, y.threshold=.005, r2.threshold=.75, tolerance=0.01)
{
  names.req <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
  check.names <- names.req %in% names(dat)
  if (sum(!check.names) > 0) 
      stop("The following required variables are missing from the dataset: ", 
          paste(names.req[!check.names], collapse = " "))
  if (sum(dat$haplo.freq) > 1 + tolerance) 
      stop("Sum of haplo freqs > 1; sum=", sum(dat$haplo.freq))
  if (sum(dat$haplo.freq) < 1 - tolerance) 
      stop("Sum of haplo freqs < 1; sum=", sum(dat$haplo.freq))
  if (abs(1 - sum(dat$haplo.freq)) > 0.01) 
      warning("Sum of haplo freqs is not 1; sum=", sum(dat$haplo.freq))
  if (sum((names(dat) %in% c("allele.freq1", "allele.freq2")))) {
      dat$allele.freq1 <- NULL
      dat$allele.freq2 <- NULL
  }
  # normalize haplotype freqs to sum to 1.0
  dat$haplo.freq <- dat$haplo.freq/sum(dat$haplo.freq)

  # generate allele freqs based on haplo freqs
  locus1 <- unique(dat$locus1)
  locus2 <- unique(dat$locus2)    
  dat$haplo.freq <- dat$haplo.freq/sum(dat$haplo.freq)
  by.vars1 <- list(dat$allele1)
  by.vars2 <- list(dat$allele2)
  names(by.vars1) <- c("allele1")
  names(by.vars2) <- c("allele2")
  af1 <- aggregate(dat$haplo.freq, by = by.vars1, FUN = sum)
  af2 <- aggregate(dat$haplo.freq, by = by.vars2, FUN = sum)
  names(af1)[length(names(af1))] <- "allele.freq1"
  names(af2)[length(names(af2))] <- "allele.freq2"
  mrg1 <- merge(dat, af1, by.x = c("allele1"), by.y = c("allele1"), all.x = T, all.y = F)
  mrg2 <- merge(mrg1, af2, by.x = c("allele2"), by.y = c("allele2"), all.x = T, all.y = F)
  dat <- mrg2

  # compute LD coefficients and expected HF under no LD
  dat$hf.exp0 <- dat$allele.freq1 * dat$allele.freq2
  dat$dij <- dat$haplo.freq - dat$hf.exp0
  allele1 <- dat$allele1
  allele2 <- dat$allele2    
  dij <- dat$dij
  exp.freq <- dat$hf.exp0    

  # convert allele names to no.s (character format) for plotting & create correspondence table
  orig.1 <- as.character(allele1) #convert from factor (or numeric)
  coded.1 <- orig.1
  for( i in 1:length(unique(orig.1)) )
  {
    coded.1[orig.1==unique(orig.1)[i]] <- paste(i) # changed from coded.1== to orig.1== for markers with allele names in the single digit sumerical range
  }
  corresp.table1 <- as.data.frame( cbind(unique(orig.1),unique(coded.1)) )
  names(corresp.table1) <- c(paste(locus1,".Allele",sep=""),paste(locus1,".Code",sep=""))

  # convert allele names to no.s (character format) for plotting & create correspondence table
  orig.2 <- as.character(allele2) #convert from factor (or numeric)
  coded.2 <- orig.2
  for( i in 1:length(unique(orig.2)) )
  {
    coded.2[orig.2==unique(orig.2)[i]] <- paste(i) # changed from coded.2== to orig.2== for markers with allele names in the single digit sumerical range
  }
  corresp.table2 <- as.data.frame( cbind(unique(orig.2),unique(coded.2)) )
  names(corresp.table2) <- c(paste(locus2,".Allele",sep=""),paste(locus2,".Code",sep=""))

  # data prep:
  al1 <- as.numeric(coded.1) #convert to numeric
  al2 <- as.numeric(coded.2)

  #plots with locus1 supplying focal allele
  count.1 <- 0; r2.1 <- NULL; npos.1 <- NULL; npos.gt.kmax.1 <- NULL; prop.pos.1 <- NULL
  maxdij.1 <- NULL; maxdij.expfreq.1 <- NULL; focal.1 <- NULL; selected.1 <- NULL
  for(i in unique(al1)) ## i <- 1 #for checking an individual plot
  {
    k <- 1 #initialize counter
    x <- NULL; y <- NULL; z <- NULL #z is the plotting character
    for(j in unique(al2)) 
    {
      # only add to x,y,z vectors if there is data for this haplo   
      if (sum(al1==i & al2==j) > 0) 
      {
        x[k] <- dij[al1==i & al2==j]
        y[k] <- exp.freq[al1==i & al2==j]
        z[k] <- j
        k <- k + 1
      }
    }
    n.neg <- sum(x<0)
    n.pos <- sum(x>0)
    n.zero <- sum(x==0)  
    n.neg.unique <- sum(unique(x)<0)
    if ( n.neg.unique>1 & n.pos>0 & max(y)>y.threshold)
    {  
      y.x.neg <- lm(y[x<0]~x[x<0])
      r2 <- summary(y.x.neg)$r.sq + 0 # fix for round() for cases where no coeff estimated
      if ( r2 > r2.threshold )
      { 
        plot(x,y,type="n",xlab="D",ylab="exp.freq", ylim=c(0,max(y)))
        title(main=paste(locus1, i ,"(",unique(orig.1)[as.numeric(i)],")",":", locus2))
        text(x,y,z)
        abline(v=0,lty=2,col=2)
        abline(h=0,lty=2,col=2)  
        abline(h=.005,lty=2,col=4)  
        title(sub=paste("R2 for D<0: ",round(r2,3)," (",unique(al2)[z[x==max(x)]],"=>",
              unique(orig.2)[z[x==max(x)]],")",sep=""))
              #z[x==max(x)] gives the coded value for allele with max dij 
        count.1 <- count.1 + 1
        r2.1[count.1] <- round(r2,3); npos.1[count.1] <- n.pos 
        maxdij.1[count.1] <- round(max(x),4); maxdij.expfreq.1[count.1] <- round(y[x==max(x)],4)
        npos.gt.kmax.1[count.1] <- sum(x>.5*max(x)) - 1
        prop.pos.1[count.1] <- round(n.pos/length(x),3)
        focal.1[count.1] <- unique(orig.1)[as.numeric(i)] 
        selected.1[count.1] <- unique(orig.2)[z[x==max(x)]] #takes 1st if there are multiple with x==max(x)
      }
    }
  }
  
  #plots with locus2 supplying focal allele
  count.2 <- 0; r2.2 <- NULL; npos.2 <- NULL; npos.gt.kmax.2 <- NULL; prop.pos.2 <- NULL
  maxdij.2 <- NULL; maxdij.expfreq.2 <- NULL; focal.2 <- NULL; selected.2 <- NULL
  for(i in unique(al2)) ## i <- 3 #for checking an individual plot
  {
    k <- 1 #initialize counter
    x <- NULL; y <- NULL; z <- NULL #z is the plotting character
    for(j in unique(al1)) 
    {
      # only add to x,y,z vectors if there is data for this haplo   
      if (sum(al2==i & al1==j) > 0)
      {
        x[k] <- dij[al2==i & al1==j]
        y[k] <- exp.freq[al2==i & al1==j]
        z[k] <- j
        k <- k + 1
      }
    }
    n.neg <- sum(x<0)
    n.pos <- sum(x>0)
    n.zero <- sum(x==0)  
    n.neg.unique <- sum(unique(x)<0)
    if ( n.neg.unique>1 & n.pos>0 & max(y)>y.threshold)
    {  
      y.x.neg <- lm(y[x<0]~x[x<0])
      r2 <- summary(y.x.neg)$r.sq + 0 # fix for round() for cases where no coeff estimated
      if ( r2 > r2.threshold )
      { 
        plot(x,y,type="n",xlab="D",ylab="exp.freq", ylim=c(0,max(y)))
        title(main=paste(locus2, i ,"(",unique(orig.2)[as.numeric(i)],")",":", locus1))
        text(x,y,z)
        abline(v=0,lty=2,col=2)
        abline(h=0,lty=2,col=2)         
        abline(h=.005,lty=2,col=4)  
        title(sub=paste("R2 for D<0: ",round(r2,3)," (",unique(al1)[z[x==max(x)]],"=>",
              unique(orig.1)[z[x==max(x)]],")",sep=""))
              #z[x==max(x)] gives the coded value for allele with max dij 
        count.2 <- count.2 + 1
        r2.2[count.2] <- round(r2,3); npos.2[count.2] <- n.pos 
        maxdij.2[count.2] <- round(max(x),4); maxdij.expfreq.2[count.2] <- round(y[x==max(x)],4)
        npos.gt.kmax.2[count.2] <- sum(x>.5*max(x)) - 1
        prop.pos.2[count.2] <- round(n.pos/length(x),3)
        focal.2[count.2] <- unique(orig.2)[as.numeric(i)] 
        selected.2[count.2] <- unique(orig.1)[z[x==max(x)]] #takes 1st if there are multiple with x==max(x)
      }   
    }
  }

  # Finish by printing correspondence tables for coded allele names
  print(corresp.table1)
  print(corresp.table2)

  focal.1 <- paste(locus1,":",focal.1,sep="")
  focal.2 <- paste(locus2,":",focal.2,sep="")
  hf.1 <- maxdij.1 + maxdij.expfreq.1
  hf.2 <- maxdij.2 + maxdij.expfreq.2
  fold.inc.1 <- round((hf.1-maxdij.expfreq.1)/maxdij.expfreq.1, 3)
  fold.inc.2 <- round((hf.2-maxdij.expfreq.2)/maxdij.expfreq.2, 3)
  summary.1 <- cbind(focal.1,selected.1,r2.1,maxdij.1,maxdij.expfreq.1,prop.pos.1,npos.gt.kmax.1,fold.inc.1)
  summary.2 <- cbind(focal.2,selected.2,r2.2,maxdij.2,maxdij.expfreq.2,prop.pos.2,npos.gt.kmax.2,fold.inc.2)
  if (dim(summary.1)[2] != 8) summary.1 <- cbind(NA,NA,NA,NA,NA,NA,NA,NA) 
  if (dim(summary.2)[2] != 8) summary.2 <- cbind(NA,NA,NA,NA,NA,NA,NA,NA) 
  summary.1 <- as.data.frame(summary.1)
  summary.2 <- as.data.frame(summary.2)
  names(summary.1) <- c("focal", "select", "r2.lt0", "maxdij", "exp.frq.max.d", "prop.gt0", "n.gt.halfmax.d", "fold.inc")
  names(summary.2) <- c("focal", "select", "r2.lt0", "maxdij", "exp.frq.max.d", "prop.gt0", "n.gt.halfmax.d", "fold.inc")
  summary.overall <- rbind(summary.1,summary.2)
  
  return(summary.overall)
}
