#!/bin/env Rscript
require(VGAM)
require(optimx)

# Read the data
data = read.table("phASER_WASP_GTEx_v8_matrix.gw_phased.txt",h=T)

# Adjustable threshold for deciding on strong ASE (1/3 is equivalent to 2-fold difference between alleles)
threshold1 = 1/3

# Adjustable threshold for deciding on moderate ASE (0.4 is equivalent to 50% greater expression from one allele)
threshold2 = 0.4

# The four components are: 1 binomial and three beta distributions (each parameterized with it's own alpha,beta values)
components = 4

write(c("gene", "area_strong", "area_moderate", "area_full"), sep=" ", 
      file = "ase_results.txt", append = TRUE, ncol=4)

# vector of initial parameter values. The parameters are (in order): the weight of the binomial (i.e. no ASE) component; 
# the weight of the first beta distribution, the weight of the second beta, alpha1, alpha2, alpha3, beta1, beta2, beta3). 
# Note: only 3 weights are given because the weights have to sum to 1 (the weight of the last component is 1 minus the sum 
# of the others)
initialvec = c(0.6, 0.37, 0.015,20.00,10.00, 2.00,20.00, 2.00,10.00)

# Adjustable threshold for deciding on strong ASE (1/3 is equivalent to 2-fold difference between alleles)
threshold1 = 1/3

# Adjustable threshold for deciding on moderate ASE (0.4 is equivalent to 50% greater expression from one allele)
threshold2 = 0.4

# The four components are: 1 binomial and three beta distributions (each parameterized with it's own alpha,beta values)
components = 4

# vector of initial parameter values. The parameters are (in order): the weight of the binomial (i.e. no ASE) component; 
# the weight of the first beta distribution, the weight of the second beta, alpha1, alpha2, alpha3, beta1, beta2, beta3). 
# Note: only 3 weights are given because the weights have to sum to 1 (the weight of the last component is 1 minus the sum 
# of the others)
#initialvec = c(0.8,0.17,0.00015,20.00,10.00,2.00,20.00,2.00,10.00)
initialvec = c(0.6, 0.37, 0.015,20.00,10.00, 2.00,20.00, 2.00,10.00)

# likelihood function
likfun <- function (pars,data,components) {
  weight = pars[1:(components-1)]
  weight = c(weight,1-sum(weight))
  alpha = pars[components:(2*components-2)]
  beta = pars[(2*components-1):length(pars)]
  if(min(weight) < 0 | min(alpha) < 0 | min(beta) < 0) { return(NA) }
  pvec = rep(0,nrow(data))
  x = data[,3] 	# counts of reads mapping to 'A' allele (If these are counts of alleles mapped to a single SNP 
  # then often the 'A' allele is the reference; sometimes the designation of which allele is 'A' and 
  # which is 'B' may be arbitrary; in other cases they may refer to specific haplotypes. This can 
  # affect the interpretation of the results.)
  n = data[,3] + data[,4] # counts of 'A' + 'B' (i.e. n is the total number of mapped reads)
  refr = 0.5 # An adjustment can be made for mapping bias in favour of the reference allele (setting to 0.5 means no adjustment) 
  pvec = weight[1] * dbinom(x,n,refr)
  for(j in 2:components) {
    pvec =  pvec + weight[j] * dbetabinom.ab(x,n,alpha[j-1],beta[j-1])
  }
  L = sum(log(pvec))
  return(-L)
}


row = 1
for (r in 1:nrow(tissue_df)){
  
  num_of_rows = as.integer(ncol(tissue_df)-1)
  
  g = data.frame(
    sample = c(1:num_of_rows),
    hap_counts = rep(0, num_of_rows),
    hap0 = rep(0, num_of_rows),
    hap1 = rep(0, num_of_rows)
  )
  
  for (x in 1:nrow(g)){
    g[as.integer(x),2] <- tissue_df[row,as.integer(x+1)]
  }
  
  colnames(g)[1] <- tissue_df[,1][row]
  
  
  n = 0
  m = 0
  k = 1 #row
  for (j in g[c(1:num_of_rows),2]){
    n = as.numeric(strsplit(j, '|', fixed=TRUE)[[1]][1])
    m = as.numeric(strsplit(j, '|', fixed=TRUE)[[1]][2])
    g[k, 3] <- as.integer(paste0(n))
    g[k, 4] <- as.integer(paste0(m))
    k = k + 1
  }
  
  
  # an adjustable parameter to optim
  max_optim_iterations = 50000
  
  # Optimize the likelihood function.
  opt = optimx(initialvec, method = "Nelder-Mead", likfun, data=g, components=components,
               control=list(maxit=max_optim_iterations))
  
  #opt = optim(initialvec,likfun,data=g,components=components, control=list(maxit=max_optim_iterations))
  
  opt_pars = c(opt$p1, opt$p2, opt$p3, opt$p4, opt$p5, opt$p6, opt$p7, opt$p8, opt$p9)
  #print(opt_pars) 	# The optimal parameter values. The first (components - 1) of these are weights. The 
  # remainder are alpha and beta parameters of the beta-binomials in your mixture distribution.
  
  pars = opt_pars
  weight = pars[1:(components - 1)]
  weight = c(weight, 1 - sum(weight))
  alpha = pars[components:(2 * components - 2)]
  beta = pars[(2 * components -1):length(pars)]
  
  
  # Next we are going to approximate the area under the probability density function. This is going to require summing over the mixture distributions
  dx = 0.001
  x = seq(dx/2, 1-dx/2, dx)
  dsum = rep(0,length(x))
  for(i in 2:components) { # summing from 2 to number of components because component 1 is the binomial (i.e. no ASE) component
    dsum = dsum + weight[i]*dbeta(x,alpha[i-1],beta[i-1])
  }
  dsum[dsum==Inf] = NA
  area_strong = 0 # area corresponding to strong ASE (i.e. this should be the proportion of the samples that show strong (> 2fold) ASE
  area_moderate = 0 # area corresponding to moderate ASE
  area_full = 0 # total area corresponding to ASE (i.e. anything that's not in the binomial (no ASE) component)
  
  for(i in 1:(length(x)-1)) {
    if(!is.na(dsum[i]) & !is.na(dsum[i+1])) {
      if((x[i] < threshold1 | x[i] > 0.5 + (0.5-threshold1))) {
        area_strong = area_strong + dx* (dsum[i]+dsum[i+1])/2
      }
      if((x[i] < threshold2 | x[i] > 0.5 + (0.5-threshold2))) {
        area_moderate = area_moderate + dx* (dsum[i]+dsum[i+1])/2
      }
      area_full = area_full + dx* (dsum[i]+dsum[i+1])/2
    }
  }
  
  row = row + 1
  
  write(c(colnames(g[1]), area_strong, area_moderate, area_full), sep=" ",
        file=paste0("ase_results.txt"), append=TRUE, ncol=4)
}
