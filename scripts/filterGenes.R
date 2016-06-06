# counts: data frame of counts by gene/sample
# threshold: level of correlation to filter out genes
excludeExtraneousGenes = function(counts, threshold = 0.96) {
  cor = cor(t(counts))
  cor[lower.tri(cor)] = NA # remove duplicates (x-y and y-x)
  cor_genes = data.frame(cbind(which(!is.na(cor), arr.ind = TRUE), na.omit(as.vector(cor))))
  cor_genes$row = rownames(cor)[cor_genes$row]
  cor_genes$col = colnames(cor)[cor_genes$col]
  cor_genes = cor_genes[cor_genes$V3 != 1 & cor_genes$V3 > threshold, ]
  if (nrow(cor_genes) == 0) return(rownames(cor))
  extraneous_genes = c()
  ok_genes = c()
  for (i in 1:nrow(cor_genes)) {
    gene1 = cor_genes$row[i]
    gene2 = cor_genes$col[i]
    if (!(gene1 %in% extraneous_genes) & !(gene1 %in% ok_genes)) {
      if (!(gene2 %in% extraneous_genes) & !(gene2 %in% ok_genes)) {
        ok_genes = c(ok_genes, gene1)
        extraneous_genes = c(extraneous_genes, gene2)
      }
      else {
        extraneous_genes = c(extraneous_genes, gene1)
      }
    }
    else {
      extraneous_genes = c(extraneous_genes, gene2)
    }
  }
  ok_genes
}

# counts: data frame of counts by gene/sample
# var: per-gene variance to keep genes
filterVariance = function(counts, var = 0.5) {
  gene_var = sapply(data.frame(t(counts)), var)
  names(gene_var) = rownames(counts)
  counts[which(gene_var > quantile(gene_var, 0.5)), ]
}
