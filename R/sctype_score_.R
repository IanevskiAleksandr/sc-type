# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs - list of gene sets positively expressed in the cell type 
# @params: gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
 
  es.max
}
