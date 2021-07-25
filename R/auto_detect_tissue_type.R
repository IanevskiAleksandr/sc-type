# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# auto_detect_tissue_type: automatically detect a tissue type of the dataset
#
# @params: path_to_db_file - DB file with cell types
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)

auto_detect_tissue_type <- function(path_to_db_file, scRNAseqData, scaled, ...){
  
  # get all tissue types in DB
  db_read = openxlsx::read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()
  
  for(tissue in tissues_){ print(paste0("Checking...", tissue));
    
    # prepare gene sets
    gs_list = gene_sets_prepare(path_to_db_file, tissue);
    
    cL_resutls = sctype_score(scRNAseqData = scRNAseqData, scaled = scaled, 
                              gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                              marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
    dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
    
    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }
  
  # order by mean score
  result_ = result_[order(-result_$score),]
  
  # plot 
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")
  
  result_
}
