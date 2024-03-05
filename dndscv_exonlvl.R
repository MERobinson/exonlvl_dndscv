#' runexondnds
#' 
#' Function to run dnds at exon level
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @author Modified by Mark Robinson for exon-level analysis
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param mutations table of mutation info
#' @param refcds exon level refcds object
#' @param outfile List of dataframes with dnds results
#' 
#' @export

runexondnds = function(mutations, refcds) {
  
  # load ref objects
  data("submod_192r_3w", package="dndscv")
  data("cancergenes_cgc81", package = "dndscv")
  refcds <- refcds[!sapply(refcds, is.null)]
  
  # set a 64x64 matrix with the functional impact of each codon transition (1=Synonymous, 2=Missense, 3=Nonsense)
  nt = c("A","C","G","T")
  trinuc_list = paste(rep(nt,each=16,times=1), rep(nt,each=4,times=4), rep(nt,each=1,times=16), sep="")
  trinuc_ind = structure(1:64, names=trinuc_list)
  trinuc_subs = NULL; for (j in 1:length(trinuc_list)) { 
    trinuc_subs = c(trinuc_subs, paste(trinuc_list[j], paste(substr(trinuc_list[j],1,1), setdiff(nt,substr(trinuc_list[j],2,2)), substr(trinuc_list[j],3,3), sep=""), sep=">")) }
  trinuc_subsind = structure(1:192, names=trinuc_subs)
  impact_matrix = array(NA, dim=c(64,64))
  colnames(impact_matrix) = rownames(impact_matrix) = trinuc_list
  for (j in 1:64) {
    for (h in 1:64) {
      from_aa = seqinr::translate(strsplit(trinuc_list[j],"")[[1]], numcode = 1)
      to_aa = seqinr::translate(strsplit(trinuc_list[h],"")[[1]], numcode = 1)
      # Annotating the impact of the mutation
      if (to_aa == from_aa){ 
        impact_matrix[j,h] = 1
      } else if (to_aa == "*"){
        impact_matrix[j,h] = 3
      } else if ((to_aa != "*") & (from_aa != "*") & (to_aa != from_aa)){
        impact_matrix[j,h] = 2
      } else if (from_aa=="*") {
        impact_matrix[j,h] = NA
      }
    }
  }
  
  ## 1. Load mutations
  message("[1/4] Parsing mutations...")
  mutations = mutations[,1:5] # Restricting input matrix to first 5 columns
  mutations[,c(1,2,3,4,5)] = lapply(mutations[,c(1,2,3,4,5)], as.character) # Factors to character
  mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
  mutations = mutations[mutations[,4]!=mutations[,5],] # Removing mutations with identical reference and mutant base
  colnames(mutations) = c("sampleID","chr","start","ref","mut")
  indna = which(is.na(mutations),arr.ind=T)
  if (nrow(indna)>0) {
    mutations = mutations[-unique(indna[,1]),] # Removing entries with an NA in any row
    warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(indna[,1]))))
  }
  mutations$end <- mutations$start + nchar(mutations$ref)
  mutations <- GRanges(mutations)
  
  ## 2. Load CDS GRanges
  message("[2/4] Preparing CDS...")
  genesgr <- lapply(refcds, "[[", "intervals_cds") %>% GRangesList()
  ol <- findOverlaps(mutations, genesgr)
  genesgr <- genesgr[unique(subjectHits(ol))]
  refcds <- refcds[unique(subjectHits(ol))]
  
  ## 3. Annotating mutations
  message("[3/4] Mapping mutational context...")
  exonref <- lapply(refcds, function(generef) {
    
    # prefilter exons
    genename <- unique(generef$gene_name)
    exonranges <- generef$intervals_cds
    ol <- findOverlaps(mutations, exonranges)
    mut <- mutations[queryHits(ol), ]
    
    # map
    gseq1up <- generef$gseq1up
    gseq1dn <- generef$gseq1dn
    cdsseq <- generef$cdsseq
    offsets <- start(mut) - start(generef$intervals_cds)[subjectHits(ol)]
    cdspos = sapply(subjectHits(ol), function(x) {
      ifelse(x == 1, 1, (sum(lengths(generef$intervals_cds)[1:(x-1)])+1)) })
    cdspos <- cdspos + offsets
    mut$cdspos <- cdspos
    mut$exonidx <- subjectHits(ol)
    mut$refanno <- sapply(cdspos, function(x) as.character(cdsseq[x]))
    mut <- mut[which(mut$ref %in% nt & mut$mut %in% nt),]
    mut <- mut[mut$ref == mut$refanno,]
    mut <- mut[mut$cdspos > 0,]
    if (length(mut) < 1) return(NULL)
    old_trinuc = sapply(seq_along(mut), function(mut_idx) {
      paste0(gseq1up[mut$cdspos][mut_idx],
             cdsseq[mut$cdspos][mut_idx],
             gseq1dn[mut$cdspos][mut_idx]) })
    new_trinuc = sapply(seq_along(mut), function(mut_idx) {
      paste0(gseq1up[mut$cdspos][mut_idx],
             mut$mut[mut_idx],
             gseq1dn[mut$cdspos][mut_idx]) })
    old_codon = as.character(Biostrings::codons(generef$cdsseq)[ceiling(mut$cdspos/3)])
    new_codon = sapply(seq_along(old_codon), function(x) {
      new = old_codon[[x]]
      substr(new, 3-mut$cdspos%%3, 3-mut$cdspos%%3) <- mut$mut[[x]]
      return(new) })
    mut_imp = impact_matrix[(trinuc_ind[new_codon]-1)*64 + trinuc_ind[old_codon]]
    mut_matrind = trinuc_subsind[paste0(old_trinuc, ">", new_trinuc)]
    
    # N for each exon
    exonref <- lapply(seq_along(generef$intervals_cds), function(exonidx) {
      if (!any(mut$exonidx == exonidx)) return(NULL)
      N = array(0, dim=c(192, 3, 1))
      idxstart <- ifelse(exonidx == 1, 1, (sum(lengths(generef$intervals_cds)[1:(exonidx-1)])*3)+1)
      idxstop <- sum(lengths(generef$intervals_cds)[1:exonidx])*3
      mut_mtrind_sub <- mut_matrind[which(mut$exonidx == exonidx)]
      mut_imp_sub <- mut_imp[which(mut$exonidx == exonidx)]
      # synonymous
      matrix_ind = table(mut_mtrind_sub[which(mut_imp_sub==1)])
      N[as.numeric(names(matrix_ind)), 1, 1] = matrix_ind
      # missense
      matrix_ind = table(mut_mtrind_sub[which(mut_imp_sub==2)])
      N[as.numeric(names(matrix_ind)), 2, 1] = matrix_ind
      # nonsense
      matrix_ind = table(mut_mtrind_sub[which(mut_imp_sub==3)])
      N[as.numeric(names(matrix_ind)), 3, 1] = matrix_ind
      list(gene_name = generef$gene_name,
           exon_id = generef$intervals_cds$exon_id[exonidx],
           exon_name = paste0(generef$gene_name, "_ex", generef$intervals_cds$exon_rank[exonidx]),
           L = generef$L[,,exonidx, drop = F],
           N = N)
    })
    return(exonref[sapply(exonref, function(x) !is.null(x))])
  })
  exonref <- flatten(exonref)
  names(exonref) <- sapply(exonref, "[[", "exon_name")
  
  ## 4. Global model
  message("[3/4] Estimating global rates...")
  
  Lall = array(sapply(exonref, function(x) x$L), dim=c(192,4,length(exonref)))
  Nall = array(sapply(exonref, function(x) x$N), dim=c(192,4,length(exonref)))
  L = apply(Lall, c(1,2), sum)
  N = apply(Nall, c(1,2), sum)
  
  # Subfunction: fitting substitution model
  fit_substmodel = function(N, L, substmodel) {
    l = c(L); n = c(N); r = c(substmodel)
    n = n[l!=0]; r = r[l!=0]; l = l[l!=0]
    params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
    indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
    colnames(indmat) = params
    for (j in 1:length(r)) {
      indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
    }
    model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
    mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
    ci = exp(confint.default(model)) # Wald confidence intervals
    par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
    return(list(par=par, model=model))
  }
  
  # Fitting all mutation rates and the 3 global selection parameters
  poissout = fit_substmodel(N, L, substmodel) # Original substitution model
  par = poissout$par
  poissmodel = poissout$model
  parmle =  setNames(par[,2], par[,1])
  mle_submodel = par
  rownames(mle_submodel) = NULL
  
  # Fitting models with 1 and 2 global selection parameters
  s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
  par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
  s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
  par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
  globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
  sel_cv = NULL
  
  # set genemut tables
  exonmuts = data.frame(exon_name = sapply(exonref, function(x) x$exon_name),
                        exon_id = sapply(exonref, function(x) x$exon_id),
                        gene_name = sapply(exonref, function(x) x$gene_name),
                        exon_number = sapply(exonref, function(x) sub("^.+_ex([0-9]+)$", "\\1", x$exon_name)),
                        n_syn=NA, n_mis=NA, n_non=NA, 
                        exp_syn=NA, exp_mis=NA, exp_non=NA, stringsAsFactors=F)
  exonmuts[,5:7] = t(sapply(exonref, function(x) colSums(x$N)))
  mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
  exonmuts[,8:10] = t(sapply(exonref, function(x) colSums(x$L*mutrates)))
  numrates = length(mutrates)
  
  ## 5. dNdScv: Negative binomial regression (with or without covariates) + local synonymous mutations
  message("[5] Running dNdScv...")
  nbreg = nbregind = NULL
  nbrdf = exonmuts[,c("n_syn","exp_syn")]
  model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
  message(sprintf("  Regression model for substitutions:theta = %0.3g.", model$theta))
  if (all(model$y==exonmuts$n_syn)) exonmuts$exp_syn_cv = model$fitted.values
  theta = model$theta
  nbreg = model
  
  # Subfunction: Analytical opt_t using only neutral subs
  mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
    tml = (n_neutral+shape-1)/(exp_rel_neutral+(1/scale))
    if (shape<=1) { # i.e. when theta<=1
      tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1, since otherwise it takes meaningless values
    }
    return(tml)
  }
  
  # Subfunction: dNdScv per gene
  selfun_cv = function(j) {
    y = as.numeric(exonmuts[j, -c(1:4)])
    x = exonref[[j]]
    exp_rel = y[4:6]/y[4]
    # Gamma
    shape = theta
    scale = y[7]/theta
    
    # a. Neutral model
    indneut = 1:3 # vector of neutral mutation types under this model (1=synonymous, 2=missense, 3=nonsense)
    opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
    mrfold = max(1e-10, opt_t/y[4]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    tmparray <- array(1, dim=c(192, 3, 1))
    ll0 = sum(dpois(x = x$N, lambda = x$L * mutrates * mrfold * tmparray, log=T)) + 
      dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik null model
    
    # b. Missense model: wmis==1, free wnon
    indneut = 1:2
    opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
    mrfold = max(1e-10, opt_t/sum(y[4])) 
    wfree = y[3]/y[6]/mrfold; wfree[y[3]==0] = 0
    tmparray <- array(1, dim=c(192, 3, 1))
    tmparray[,3,] <- wfree
    llmis = sum(dpois(x = x$N, lambda = x$L * mutrates * mrfold * tmparray, log=T)) + 
      dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
    
    # c. Nonsense muts model: free wmis, wnon==1
    indneut = c(1,3)
    opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
    mrfold = max(1e-10, opt_t/sum(y[4])) 
    wfree = y[2]/y[5]/mrfold; wfree[y[2]==0] = 0
    tmparray <- array(1, dim=c(192, 3, 1))
    tmparray[,2,] <- wfree
    llnon = sum(dpois(x = x$N, lambda = x$L * mutrates * mrfold * tmparray, log=T)) + 
      dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wnon
    
    # d. Free selection model: free wmis, free wnon
    indneut = 1
    opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
    mrfold = max(1e-10, opt_t/sum(y[4])) 
    wfree = y[2:3]/y[5:6]/mrfold; wfree[y[2:3]==0] = 0
    tmparray <- array(1, dim=c(192, 3, 1))
    tmparray[,2:3,] <- wfree
    llall_unc = sum(dpois(x = x$N, lambda = x$L * mutrates * mrfold * tmparray, log=T)) + 
      dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free both
    tmparray <- array(1, dim=c(192, 3, 1))
    wmisfree = y[2]/y[5]/mrfold; wmisfree[y[2]==0] = 0
    tmparray[,2,] <- wmisfree
    wnonfree = sum(y[3])/sum(y[6])/mrfold; wnonfree[sum(y[3])==0] = 0
    tmparray[,3,] <- wnonfree
    llall = sum(dpois(x = x$N, lambda= x$L * mutrates * mrfold * tmparray, log=T)) + 
      dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis, free wnon
    p = 1-pchisq(2*c(llall_unc-llmis, llall-c(llnon,ll0)), df=c(1,1,2))
    return(c(wmisfree, wnonfree, p))
  }
  sel_cv = as.data.frame(t(sapply(1:nrow(exonmuts), selfun_cv)))
  colnames(sel_cv) = c("wmis_cv","wnon_cv","pmis_cv","pnon_cv","pallsubs_cv")
  sel_cv$qmis_cv = p.adjust(sel_cv$pmis_cv, method="BH")
  sel_cv$qnon_cv = p.adjust(sel_cv$pnon_cv, method="BH")
  sel_cv$qallsubs_cv = p.adjust(sel_cv$pallsubs_cv, method="BH")
  sel_cv = cbind(exonmuts[,1:7], sel_cv)
  sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$pnon_cv, -sel_cv$wmis_cv),] 
  list(globaldnds = globaldnds, 
       sel_cv = sel_cv)
}

