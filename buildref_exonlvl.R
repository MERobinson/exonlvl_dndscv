#' builexondref
#' 
#' Function to build a RefCDS object from a reference genome and a table of transcripts. 
#' The RefCDS object has to be precomputed for any new species or assembly prior to running dndscv. 
#' This function generates an .rda file that needs to be input into dndscv using the refdb argument. 
#' Note that when multiple CDS share the same gene name (second column of cdsfile), 
#' the longest coding CDS will be chosen for the gene. CDS with ambiguous bases (N) will not be considered.
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @author Modified by Mark Robinson for exon-level analysis
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param txdb EnsDB transcript database
#' @param bsgenome BSgenome object for genome build matching EnsDB
#' @param outfile Output file name (default = "RefCDS.rda").
#' 
#' @export

builexondref = function(txdb, bsgenome, outfile = "RefCDS.rda") {
  
  ## 1. Load & filter CDS GRanges
  message("[1/2] Preparing CDS...")
  genetx <- transcriptsBy(txdb, "gene", filter = TxBiotypeFilter("protein_coding"))
  seqlevelsStyle(genetx) <- "UCSC"
  genome(genetx) <- "hg38"
  isCircular(genetx)["chrM"] <- TRUE
  genetx <- updateObject(keepStandardChromosomes(genetx, pruning.mode = "coarse"))
  cdsgr <- cdsBy(txdb, "tx", filter = TxBiotypeFilter("protein_coding"))
  seqlevelsStyle(cdsgr) <- "UCSC"
  genome(cdsgr) <- "hg38"
  isCircular(cdsgr)["chrM"] <- TRUE
  cdsgr <- updateObject(keepStandardChromosomes(cdsgr, pruning.mode = "coarse"))
  
  # Precalculating a 64x64 matrix with the functional impact of each codon transition (1=Synonymous, 2=Missense, 3=Nonsense)
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
  
  ## 2. Build the RefCDS object
  message("[2/2] Building the RefCDS object...")
  RefCDS <- lapply(genetx, function(tx) {
    
    # get CDS seq & flanking bases
    geneid <- unique(tx$gene_id)
    txcds <- cdsgr[names(cdsgr) %in% tx$tx_name]
    cdslengths <- sapply(txcds, function(y) sum(lengths(y)))
    txcds <- txcds[names(cdslengths)[order(-cdslengths)],]
    generef <- lapply(txcds, function(cds) {
      strand <- as.character(strand(cds))[1]
      chr <- as.character(seqnames(cds))[1]
      cdsseq <- getSeq(bsgenome, cds+1)
      gseq1dn <- sapply(cdsseq, "[", -c(1:2)) %>% DNAStringSet()
      gseq1up <- sapply(cdsseq, function(x) x[-c((length(x)-1):(length(x)))]) %>% DNAStringSet()
      cdsseq <- sapply(cdsseq, function(x) x[-c(1,length(x))]) %>% DNAStringSet()
      gseq1dn <- DNAString(paste(as.character(gseq1dn), collapse = ""))
      gseq1up <- DNAString(paste(as.character(gseq1up), collapse = ""))
      cdsseq <- DNAString(paste(as.character(cdsseq), collapse = ""))
      pseq <- suppressWarnings(translate(cdsseq))
      if (grepl('*', pseq[-length(pseq)], fixed = T) | grepl("N", cdsseq) | length(cdsseq)%%3 != 0) return(NULL) 
      list(gene_name = unname(mapIds(txdb, keys = geneid, column = "SYMBOL", keytype = "GENEID")),
           CDS_length = length(cdsseq),
           chr = chr,
           strand = strand,
           intervals_cds = cds,
           cdsseq = cdsseq,
           gseq1up = gseq1up,
           gseq1dn = gseq1dn)
    })
    generef <- generef[sapply(generef, function(x) !is.null(x))]
    if (length(generef) < 1) { print(paste0("No valid CDS found for gene ", geneid)); return(NULL) }
    generef <- generef[[1]]
    
    # set separate L for each exon
    cdsseq <- as.character(as.vector(generef$cdsseq))
    gseq1up <- as.character(as.vector(generef$gseq1up))
    gseq1dn <- as.character(as.vector(generef$gseq1dn))
    ind = rep(1:length(cdsseq), each=3)
    old_trinuc = paste(gseq1up[ind], cdsseq[ind], gseq1dn[ind], sep="")
    new_base = c(sapply(cdsseq, function(x) nt[nt!=x]))
    new_trinuc = paste(gseq1up[ind], new_base, gseq1dn[ind], sep="")
    codon_start = rep(seq(1,length(cdsseq),by=3),each=9)
    old_codon = paste(cdsseq[codon_start], cdsseq[codon_start+1], cdsseq[codon_start+2], sep="")
    pos_in_codon = rep(rep(1:3, each=3), length.out=length(old_codon))
    aux = strsplit(old_codon,"")
    new_codon = sapply(1:length(old_codon), function(x) { new_codonx = aux[[x]]; new_codonx[pos_in_codon[x]] = new_base[x]; return(new_codonx) } )
    new_codon = paste(new_codon[1,], new_codon[2,], new_codon[3,], sep="")
    imp = impact_matrix[(trinuc_ind[new_codon]-1)*64 + trinuc_ind[old_codon]]
    matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
    L = array(0, dim=c(192, 3, length(generef$intervals)))
    for (exonidx in 1:length(generef$intervals)) {
      idxstart <- ifelse(exonidx == 1, 1, (sum(lengths(generef$intervals)[1:(exonidx-1)])*3)+1)
      idxstop <- sum(lengths(generef$intervals)[1:exonidx])*3
      mtrindsub <- matrind[idxstart:idxstop]
      impsub <- imp[idxstart:idxstop]
      # synonymous
      matrix_ind = table(mtrindsub[which(impsub==1)])
      L[as.numeric(names(matrix_ind)), 1, exonidx] = matrix_ind
      # missense
      matrix_ind = table(mtrindsub[which(impsub==2)])
      L[as.numeric(names(matrix_ind)), 2, exonidx] = matrix_ind
      # nonsense
      matrix_ind = table(mtrindsub[which(impsub==3)])
      L[as.numeric(names(matrix_ind)), 3, exonidx] = matrix_ind
    }
    generef$L <- L
    return(generef)
  })
  save(RefCDS, file=outfile)
  return(RefCDS)
}
