#' Constrained Disequilibrium Values (CDV).
#'
#' A function to perform a Constrained Disequilibrium Values (CDV) analysis for 3-locus haplotype frequency data. 
#' 
#' @param dat	A data.frame with 5 required variables (having the names listed below):
#'      \tabular{ll}{
#'           \code{haplo.freq} \tab A numeric vector of haplotype frequencies.\cr
#'            \code{locus1} \tab A character vector indentifying the first locus.\cr
#'            \code{locus2} \tab A character vector indentifying the second locus.\cr
#'            \code{allele1} \tab A character vector indentifying the allele at locus 1.\cr
#'            \code{allele2} \tab A character vector indentifying the allele at locus 2.\cr
#'            }
#'      if the variables \code{allele.freq1} and \code{allele.freq2} are on the data file, the
#'      it is assumed that the data may contain an incomplete set of haplotypes whose frequencies
#'      do not sum to 1.0 (and thus allele freqs can not be computed from haplotype freqs). The
#'      supplied allele.freqs are used in the CDV computations instead.
#'
#' @param tolerance A threshold for the sum of the haplotype frequencies.
#'    If the sum of the haplotype frequencies is greater than 1+tolerance or less
#'    than 1-tolerance an error is returned. The default is 0.01.
#'
#' @return The return value is a dataframe with the following components:
#'  \tabular{ll}{
#'  \code{trio} 	\tab set of 3 alleles\cr
#'  \code{l0}   	\tab the conditioned upon or constraining allele\cr
#'  \code{hf}   	\tab the haplotype freq for the haplo consisting of the other two loci\cr
#'  \code{d}      \tab the ordinary disequilibrium coefficient \cr
#'  \code{dprime}	\tab the normalized disequilibrium coefficient (d') \cr
#'  \code{dprime2}\tab the normalized diseq coeff (d'') considering constraints by third locus\cr
#'  \code{delta}	\tab |d' - d''| \cr
#'  }
#'
#' @section Details:
#' A warning message is given if the sum of the haplotype frequencies is greater than 1.01 or less
#' than 0.99 (regardless of the \code{tolerance} setting). The haplotype frequencies that are passed
#' to the function are normalized within the function to sum to 1.0 by dividing each frequency by 
#' the sum of the passed frequencies.
#'
#' @examples
#' # An example using haplotype frequencies from Northern Ireland
#' library(LDtools)
#' data(NIreland.freqs)
#' ni.dat <- NIreland.freqs
#' abd.dat  = ni.dat[ni.dat$locus1 %in% c("A","B","DRB1") & ni.dat$locus2 %in% c("A","B","DRB1"),]
#' abd.cdv = compute.CDV(abd.dat)
#' abd.out = flag.CDVselect(abd.cdv)
#' abd.out[abd.out$flag.CDV >0 & abd.out$flag.hf==1,]
#' 
#' # Note that three haplotypes are identified as having experienced selection based on 
#' # D' (dprime) and D'' (dprime2) values: A*01:01~B*08:01~DRB1*03:01 (with HLA-B*08:01 as
#' # the putatively selected allele); A*29:02~B*44:03~DRB1*07:01 (with HLA-B*44:03 as the 
#' # putatively selected allele); and A*31:01~B*40:01~DRB1*04:04 (with HLA-A*31:01 as the 
#' # putatively selected allele)
#' @export
#' @importFrom stats aggregate

compute.LDfromHF <- function(dat, tolerance=0.01)
{
  loci <- unique(c(dat$locus1,dat$locus2))
  if (sum(dat$haplo.freq) > 1 + tolerance) 
      stop("Sum of haplo freqs > 1; sum=", sum(dat$haplo.freq), " for loci ", loci)

  if (sum(dat$haplo.freq) < 1 - tolerance) 
      stop("Sum of haplo freqs < 1; sum=", sum(dat$haplo.freq), " for loci ", loci)
      
  if (abs(1 - sum(dat$haplo.freq)) > 0.01) 
      warning("Sum of haplo freqs is not 1; sum=", sum(dat$haplo.freq), " for loci", loci)
  
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
  mrg1 <- merge(dat, af1, by.x = c("allele1"), by.y = c("allele1"), all.x = TRUE, all.y = FALSE)
  mrg2 <- merge(mrg1, af2, by.x = c("allele2"), by.y = c("allele2"), all.x = TRUE, all.y = FALSE)
  dat <- mrg2

  # compute LD coefficients
  dat$d <- dat$haplo.freq - dat$allele.freq1 * dat$allele.freq2
  dat$dprime.den <- rep(0,dim(dat)[1])
  dat$den.lt0 <- pmin( dat$allele.freq1*dat$allele.freq2, (1-dat$allele.freq1)*(1-dat$allele.freq2) )
  dat$den.ge0 <- pmin( (1-dat$allele.freq1)*dat$allele.freq2, dat$allele.freq1*(1-dat$allele.freq2) )
  dat$dprime.den[dat$d < 0] <- dat$den.lt0[dat$d < 0]
  dat$dprime.den[dat$d >=0] <- dat$den.ge0[dat$d >=0]
  dat$dprime <- dat$d/dat$dprime.den  
  dat$dprime[dat$d ==0] <- 0 ########
  dat$den.lt0 <- NULL
  dat$den.ge0 <- NULL
  dat$dprime.den <- NULL
  dat$hf.exp0 <- NULL
  
  allele1 <- dat$allele1
  allele2 <- dat$allele2    
  
  return(dat)
}

compute.CDV <- function(dat, tolerance=0.01)
{
  my.eps = .0001
  names.req <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
  check.names <- names.req %in% names(dat)
  if (sum(!check.names) > 0) 
      stop("The following required variables are missing from the dataset: ", 
          paste(names.req[!check.names], collapse = " "))
  partial.haplos = FALSE
  if (sum((names(dat) %in% c("allele.freq1", "allele.freq2")))==2)
      partial.haplos = TRUE
  if (sum((names(dat) %in% c("allele.freq1", "allele.freq2")))==1){
      dat$allele.freq1 <- NULL
      dat$allele.freq2 <- NULL
  } 
  loci <- unique(c(dat$locus1,dat$locus2))
  if (length(loci) > 3) 
      stop("There are more than 3 loci on the data file")

  if (partial.haplos)
  {
    dat$d <- dat$haplo.freq - dat$allele.freq1 * dat$allele.freq2
    dat$dprime.den <- rep(0,dim(dat)[1])
    dat$den.lt0 <- pmin( dat$allele.freq1*dat$allele.freq2, (1-dat$allele.freq1)*(1-dat$allele.freq2) )
    dat$den.ge0 <- pmin( (1-dat$allele.freq1)*dat$allele.freq2, dat$allele.freq1*(1-dat$allele.freq2) )
    dat$dprime.den[dat$d < 0] <- dat$den.lt0[dat$d < 0]
    dat$dprime.den[dat$d >=0] <- dat$den.ge0[dat$d >=0]
    dat$dprime <- dat$d/dat$dprime.den  
    dat$den.lt0 <- NULL
    dat$den.ge0 <- NULL
    dat$dprime.den <- NULL
  }
    
  dat12 <- dat[dat$locus1==loci[1] & dat$locus2==loci[2],]
  dat12b<- dat[dat$locus1==loci[2] & dat$locus2==loci[1],]
  if (dim(dat12b)[1]>0 & dim(dat12)[1]>0) stop("Problem with data for 1st & 2nd locus")
  if (dim(dat12b)[1]>0 & dim(dat12)[1]==0) dat12 <- dat12b
  dat13 <- dat[dat$locus1==loci[1] & dat$locus2==loci[3],]
  dat13b<- dat[dat$locus1==loci[3] & dat$locus2==loci[1],]
  if (dim(dat13b)[1]>0 & dim(dat13)[1]>0) stop("Problem with data for 1st & 3rd locus")
  if (dim(dat13b)[1]>0 & dim(dat13)[1]==0) dat13 <- dat13b
  dat23 <- dat[dat$locus1==loci[2] & dat$locus2==loci[3],]
  dat23b<- dat[dat$locus1==loci[3] & dat$locus2==loci[2],]
  if (dim(dat23b)[1]>0 & dim(dat23)[1]>0) stop("Problem with data for 2nd & 3rd locus")
  if (dim(dat23b)[1]>0 & dim(dat23)[1]==0) dat23 <- dat23b

  if (!partial.haplos)
  {
    dat12b <- compute.LDfromHF(dat12, tolerance)
    dat13b <- compute.LDfromHF(dat13, tolerance)
    dat23b <- compute.LDfromHF(dat23, tolerance)
    dat12 = dat12b; rm(dat12b)
    dat13 = dat13b; rm(dat13b)
    dat23 = dat23b; rm(dat23b)
  } else
  {
    if (abs(1 - sum(dat12$haplo.freq)) > 0.01) 
       warning("Assuming partial haplotypes supplied - Sum of haplo freqs =", sum(dat12$haplo.freq), " for loci ", unique(c(dat12$locus1,dat12$locus2)))
    if (abs(1 - sum(dat13$haplo.freq)) > 0.01) 
       warning("Assuming partial haplotypes supplied - Sum of haplo freqs =", sum(dat13$haplo.freq), " for loci ", unique(c(dat13$locus1,dat13$locus2)))
    if (abs(1 - sum(dat23$haplo.freq)) > 0.01) 
       warning("Assuming partial haplotypes supplied - Sum of haplo freqs =", sum(dat23$haplo.freq), " for loci ", unique(c(dat23$locus1,dat23$locus2)))
  }
 #dat <- do.call("rbind", list(dat12, dat13, dat23))
  
  out.list1 = sort( unique(c(dat12$allele1,dat13$allele1)) )
  out.list2 = sort( unique(c(dat12$allele2,dat23$allele1)) ) #Note: allele2 for dat12
  out.list3 = sort( unique(c(dat13$allele2,dat23$allele2)) ) #Note: allele2 for both

  rbind.out1 = NULL
  rbind.out2 = NULL
  rbind.out3 = NULL
  for (out.locus in loci[1:3])
  {
    if (out.locus==loci[1]) out.list = out.list1
    if (out.locus==loci[2]) out.list = out.list2
    if (out.locus==loci[3]) out.list = out.list3
    for (out.allele in out.list)
    {
      if (out.locus==loci[1])
      {
        tmp2 = dat12[dat12$allele1==out.allele,] #dat12$locus1==out.locus
        tmp3 = dat13[dat13$allele1==out.allele,] #dat13$locus1==out.locus
        p1a = unique(tmp2$allele.freq1)
        p1b = unique(tmp3$allele.freq1)
        if (length(p1a)==1 & length(p1b)==1) p1 = mean(p1a,p1b)
        if (abs(p1a-p1b)>tolerance) stop("ERROR: abs(p1a-p1b)>tolerance)")
        a2.list = sort(unique(tmp2$allele2))
        a3.list = sort(unique(tmp3$allele2))
        for (a2 in a2.list)
        {
          p2 = tmp2$allele.freq2[tmp2$allele2==a2]
          for (a3 in a3.list)
          {
            p3 = tmp3$allele.freq2[tmp3$allele2==a3]
            q1=1-p1; q2=1-p2; q3=1-p3
            d12 = tmp2$d[tmp2$allele2==a2]
            d13 = tmp3$d[tmp3$allele2==a3]
            d23 = dat23$d[dat23$allele1==a2 & dat23$allele2==a3]
            hf = dat23$haplo.freq[dat23$allele1==a2 & dat23$allele2==a3]
           #if (length(hf)==0) hf = 9
            if (length(hf)>0)
            {
              if (length(d12)==0) d12=0-p1*p2
              if (length(d13)==0) d13=0-p1*p3
              if (length(d23)==0) d23=0-p2*p3
              #out.locus=1 (c->1; a->2; b->3) 
              maxD1 = p2 * q3
              maxD2 = q2 * p3
              maxD3 = p2 * q3 * p1 + q2 * p3 * q1 + d12 - d13    #NB: d12=d21 #keep  
              maxD4 = p2 * q3 * q1 + q2 * p3 * p1 - d12 + d13    #NB: d13=d31 #keep
              minD1 = p2 * p3 * (-1.0)
              minD2 = q2 * q3 * (-1.0)
              minD3 = (p2 * p3 * p1 + q2 * q3 * q1 + d12 + d13) * (-1.0) #keep
              minD4 = (p2 * p3 * q1 + q2 * q3 * p1 - d12 - d13) * (-1.0) #keep
              maxD = min(c(maxD1,maxD2,maxD3,maxD4))
              minD = max(c(minD1,minD2,minD3,minD4)) 
              #-------------------
              #if (d23>0 & minD<=0) dprime2 = d23/maxD
              #if (d23>0 & minD>0) dprime2 = (d23-minD)/(maxD-minD)
              #if (d23<0 & maxD>=0) dprime2 = d23/(-1*minD)
              #if (d23<0 & maxD<0) dprime2 = (d23-maxD)/(maxD-minD)
              #-------------------
              if (d23>0 & minD<=0) 
              {
                if (abs(maxD)<my.eps) dprime2 = 0.0
                else dprime2 = d23/maxD
              }  
              if (d23>0 & minD>0) 
              { 
                if (abs(maxD-minD)<my.eps) dprime2 = 1.0
                else dprime2 = (d23-minD)/(maxD-minD)
              }
              if (d23<0 & maxD>=0) 
              {
                if (abs(minD)<my.eps) dprime2 = 0.0
                else dprime2 = d23/(-1*minD)
              } 
              if (d23<0 & maxD<0) 
              {
                if (abs(maxD-minD)<my.eps) dprime2 = -1.0
                else dprime2 = (d23-maxD)/(maxD-minD)
              }
              if (d23==0.0) dprime2 = 0.0
              dprime = dat23$dprime[dat23$allele1==a2 & dat23$allele2==a3]
              if (length(dprime)==0) 
              { 
                if (d23==0) dprime = 0
                if (d23> 0) dprime = d23/min(p2*q2,q2*p3)            #****
                if (d23< 0) dprime = d23/(-1*max(-1*p2*p3,-1*q2*q3)) #****
              }
              out = c(out.locus, out.allele, unique(tmp2$locus2), a2, unique(tmp3$locus2), a3, dprime, dprime2, p1, p2, p3, d23, hf)
              rbind.out1 = rbind(rbind.out1, out)
            } #if (length(hf)>0)
          } #for(a3)
        } #for(a2)
      } #if(out.locus)
      #================================================
      if (out.locus==loci[2])
      {
        tmp1 = dat12[dat12$allele2==out.allele,] #dat12$locus2==out.locus
        tmp3 = dat23[dat23$allele1==out.allele,] #dat23$locus1==out.locus
        p2a = unique(tmp1$allele.freq2)
        p2b = unique(tmp3$allele.freq1)
        if (length(p2a)==1 & length(p2b)==1) p2 = mean(p2a,p2b)
        if (abs(p2a-p2b)>tolerance) stop("ERROR: abs(p2a-p2b)>tolerance)")
        a1.list = sort(unique(tmp1$allele1))
        a3.list = sort(unique(tmp3$allele2))
        for (a1 in a1.list)
        {
          p1 = tmp1$allele.freq1[tmp1$allele1==a1]
          for (a3 in a3.list)
          {
            p3 = tmp3$allele.freq2[tmp3$allele2==a3]
            q1=1-p1; q2=1-p2; q3=1-p3
            d12 = tmp1$d[tmp1$allele1==a1]
            d13 = dat13$d[dat13$allele1==a1 & dat13$allele2==a3]
            d23 = tmp3$d[tmp3$allele2==a3]
            hf = dat13$haplo.freq[dat13$allele1==a1 & dat13$allele2==a3]
            #if (length(hf)==0) hf = 9
            if (length(hf)>0)
            {
              if (length(d12)==0) d12=0-p1*p2
              if (length(d13)==0) d13=0-p1*p3
              if (length(d23)==0) d23=0-p2*p3
              #out.locus=2 (c->2; a->1; b->3) 
              maxD1 = p1 * q3
              maxD2 = q1 * p3
              maxD3 = p1 * q3 * p2 + q1 * p3 * q2 + d12 - d23    
              maxD4 = p1 * q3 * q2 + q1 * p3 * p2 - d12 + d23    #NB: d23=d32
              minD1 = p1 * p3 * (-1.0)
              minD2 = q1 * q3 * (-1.0)
              minD3 = (p1 * p3 * p2 + q1 * q3 * q2 + d12 + d23) * (-1.0)
              minD4 = (p1 * p3 * q2 + q1 * q3 * p2 - d12 - d23) * (-1.0)
              maxD = min(c(maxD1,maxD2,maxD3,maxD4))
              minD = max(c(minD1,minD2,minD3,minD4)) 
              #-------------------
              #if (d13>0 & minD<=0) dprime2 = d13/maxD
              #if (d13>0 & minD>0) dprime2 = (d13-minD)/(maxD-minD)
              #if (d13<0 & maxD>=0) dprime2 = d13/(-1*minD)
              #if (d13<0 & maxD<0) dprime2 = (d13-maxD)/(maxD-minD)
              #-------------------
              if (d13>0 & minD<=0) 
              {
                if (abs(maxD)<my.eps) dprime2 = 0.0
                else dprime2 = d13/maxD
              }  
              if (d13>0 & minD>0) 
              { 
                if (abs(maxD-minD)<my.eps) dprime2 = 1.0
                else dprime2 = (d13-minD)/(maxD-minD)
              }
              if (d13<0 & maxD>=0) 
              {
                if (abs(minD)<my.eps) dprime2 = 0.0
                else dprime2 = d13/(-1*minD)
              } 
              if (d13<0 & maxD<0) 
              {
                if (abs(maxD-minD)<my.eps) dprime2 = -1.0
                else dprime2 = (d13-maxD)/(maxD-minD)
              }
              if (d13==0.0) dprime2 = 0.0
              dprime = dat13$dprime[dat13$allele1==a1 & dat13$allele2==a3]
              if (length(dprime)==0) 
              { 
                if (d13==0) dprime = 0
                if (d13> 0) dprime = d13/min(p1*q3,q1*p3)            #****
                if (d13< 0) dprime = d13/(-1*max(-1*p1*p3,-1*q1*q3)) #****
              }
              out = c(out.locus, out.allele, unique(tmp1$locus1), a1, unique(tmp3$locus2), a3, dprime, dprime2, p2, p1, p3, d13, hf)
              rbind.out2 = rbind(rbind.out2, out)
            } #if (length(hf)>0)
          } #for(a3)
        } #for(a1)
      } #if(out.locus)
      #================================================
      if (out.locus==loci[3])
      {
        tmp1 = dat13[dat13$allele2==out.allele,] #dat13$locus2==out.locus
        tmp2 = dat23[dat23$allele2==out.allele,] #dat23$locus2==out.locus
        p3a = unique(tmp1$allele.freq2) #NB: allele.freq2 for out.locus
        p3b = unique(tmp2$allele.freq2)
        if (length(p3a)==1 & length(p3b)==1) p3 = mean(p3a,p3b)
        if (abs(p3a-p3b)>tolerance) stop("ERROR: abs(p3a-p3b)>tolerance)")
        a1.list = sort(unique(tmp1$allele1))
        a2.list = sort(unique(tmp2$allele1))
        for (a1 in a1.list)
        {
          p1 = tmp1$allele.freq1[tmp1$allele1==a1]
          for (a2 in a2.list)
          {
            p2 = tmp2$allele.freq1[tmp2$allele1==a2]
            q1=1-p1; q2=1-p2; q3=1-p3
            d13 = tmp1$d[tmp1$allele1==a1]
            d23 = tmp2$d[tmp2$allele1==a2]
            d12 = dat12$d[dat12$allele1==a1 & dat12$allele2==a2]
            hf = dat12$haplo.freq[dat12$allele1==a1 & dat12$allele2==a2]
            #if (length(hf)==0) hf = 9
            if (length(hf)>0)
            {
              if (length(d12)==0) d12=0-p1*p2
              if (length(d13)==0) d13=0-p1*p3
              if (length(d23)==0) d23=0-p2*p3
              #out.locus=3 (c->3; a->1; b->2)
              maxD1 = p1 * q2
              maxD2 = q1 * p2
              maxD3 = p1 * q2 * p3 + q1 * p2 * q3 + d13 - d23
              maxD4 = p1 * q2 * q3 + q1 * p2 * p3 - d13 + d23
              minD1 = p1 * p2 * (-1.0)
              minD2 = q1 * q2 * (-1.0)
              minD3 = (p1 * p2 * p3 + q1 * q2 * q3 + d13 + d23) * (-1.0)
              minD4 = (p1 * p2 * q3 + q1 * q2 * p3 - d13 - d23) * (-1.0)
              maxD = min(c(maxD1,maxD2,maxD3,maxD4))
              minD = max(c(minD1,minD2,minD3,minD4)) 
              #-------------------
              #if (d12>0 & minD<=0) dprime2 = d12/maxD
              #if (d12>0 & minD>0) dprime2 = (d12-minD)/(maxD-minD)
              #if (d12<0 & maxD>=0) dprime2 = d12/(-1*minD)
              #if (d12<0 & maxD<0) dprime2 = (d12-maxD)/(maxD-minD)
              #-------------------
              if (d12>0 & minD<=0) 
              {
                if (abs(maxD)<my.eps) dprime2 = 0.0
                else dprime2 = d12/maxD
              }  
              if (d12>0 & minD>0) 
              { 
                if (abs(maxD-minD)<my.eps) dprime2 = 1.0
                else dprime2 = (d12-minD)/(maxD-minD)
              }
              if (d12<0 & maxD>=0) 
              {
                if (abs(minD)<my.eps) dprime2 = 0.0
                else dprime2 = d12/(-1*minD)
              } 
              if (d12<0 & maxD<0) 
              {
                if (abs(maxD-minD)<my.eps) dprime2 = -1.0
                else dprime2 = (d12-maxD)/(maxD-minD)
              }
              if (d12==0.0) dprime2 = 0.0
              #-------------------
              dprime = dat12$dprime[dat12$allele1==a1 & dat12$allele2==a2]
              if (length(dprime)==0) 
              { 
                if (d12==0) dprime = 0
                if (d12> 0) dprime = d12/min(p1*q1,q1*p2)
                if (d12< 0) dprime = d12/(-1*max(-1*p1*p2,-1*q1*q2))
              }
              out = c(out.locus, out.allele, unique(tmp1$locus1), a1, unique(tmp2$locus1), a2, dprime, dprime2, p3, p1, p2, d12, hf)
              rbind.out3 = rbind(rbind.out3, out)
            } #if (length(hf)>0)
          } #for(a2)
        } #for(a1)
      } #if(out.locus)
      #================================================          
    } #for(out.allele)
  } #for(out.locus)   
  rbind.out = rbind(rbind.out1, rbind.out2, rbind.out3)
  out.df = as.data.frame(rbind.out)
  rownames(out.df) = NULL
  colnames(out.df) = c("out.loc","out.a","loc1","allele1","loc2","allele2","dprime","dprime2", "af.out", "af.1", "af.2", "d", "hf")
  out.df$dprime = as.numeric(as.character(out.df$dprime))
  out.df$dprime2= as.numeric(as.character(out.df$dprime2))
  out.df$dprime[out.df$dprime >  1] =  1
  out.df$dprime[out.df$dprime < -1] = -1
  out.df$dprime2[out.df$dprime2 >  1] =  1
  out.df$dprime2[out.df$dprime2 < -1] = -1
  out.df$af.1 = as.numeric(as.character(out.df$af.1))
  out.df$af.2 = as.numeric(as.character(out.df$af.2))
  out.df$af.out = as.numeric(as.character(out.df$af.out))
  out.df$delta = abs(out.df$dprime) - abs(out.df$dprime2)
  out.df$hf = as.numeric(as.character(out.df$hf))

  out.df$l0 = paste(out.df$out.loc, out.df$out.a, sep="_")
  out.df$l1 = paste(out.df$loc1, out.df$allele1, sep="_")
  out.df$l2 = paste(out.df$loc2, out.df$allele2, sep="_")
  tmp = as.data.frame( t(apply(out.df[,c("l0","l1","l2")], 1, sort, )) )
  names(tmp) = c("n1","n2","n3")
  out.df$trio = paste(tmp$n1,tmp$n2,tmp$n3, sep=":")
  trio.n = table(out.df$trio)
  trio.n.3 = names(trio.n[trio.n==3])
  out.df = out.df[out.df$trio %in% trio.n.3, c("trio","l0","hf","d","dprime","dprime2","delta")]
  out.df = lsort(out.df, by="trio")
  return(out.df)  
} #compute.CDV  
#========================================================
flag.CDVselect <- function(out, hf.min=0.01)
{
  tol2=.0001
  n.trio = dim(out)[1]/3
  i=1; flag.delta = NULL; flag.hf = NULL
  for (trio in 1:n.trio - 1)
  { 
    flag=0
    hf.min = min(out$hf[i], out$hf[i+1], out$hf[i+2])
    if (hf.min >= .01) flag.hf = c(flag.hf, rep(1,3)) 
    else flag.hf = c(flag.hf, rep(0,3)) 
    d1 = out$delta[i]
    d2 = out$delta[i+1]
    d3 = out$delta[i+2]
    d.max = max(c(d1,d2,d3))
    d.med = median(c(d1,d2,d3))
    if (d1<=0 & d2<=0 & d3<=0) flag=0 #otherwise, >=1 is pos
    else {
      if (d1 <=0 & d2<=0 | d1 <=0 & d3<=0 | d2 <=0 & d3<=0) {
        if (d1<=0 & d2<=0 & d.max>tol2) flag=1.3
        if (d1<=0 & d3<=0 & d.max>tol2) flag=1.2 
        if (d2<=0 & d3<=0 & d.max>tol2) flag=1.1
      } else if (d.max > 2*d.med ) {
        if (d1 == d.max & d.max>tol2) flag=2.1
        if (d2 == d.max & d.max>tol2) flag=2.2
        if (d3 == d.max & d.max>tol2) flag=2.3
      } 
    }
    flag.delta = c(flag.delta, rep(flag,3))
    i = i + 3
  }
  out$flag.CDV = flag.delta
  out$flag.hf = flag.hf
  return(out)  
}  

