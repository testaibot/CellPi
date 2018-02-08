


library(scater)     #main data container
library(scran)      #zero-inflated normalization by pooling
library(dropestr)   #UMI correction
library(SAVER)      #impute dropouts step
library(doParallel) #required by SAVER
library(monocle)    #census relative2abs used for "pseudo-spike-ins" normalization
library(biomaRt)    #genome annotations
library(M3Drop)     #find DE genes
library(Seurat)     #unsupervised clustering

library(VennDiagram)
library(gplots)

library(clusterProfiler)#GO analysis stage annotation convertercompatible with ReactomePA
library(ReactomePA)#GO and GSEA





types2names = function(types)
{
  dict = rep(1, length(unique(types)))
  names(dict) = unique(types)
  
  cell_names = NULL
  for (type in types)
  {
    cell_names = c(cell_names, paste0(type, ".", dict[type]))
    dict[type] = dict[type] + 1
  }
  return(cell_names)
}

annot_path = function(spec)
{
  if(spec %in% c("h", "human"))
  {
    return("human_cycle.rds")
  }
  if(spec %in% c("m", "mouse"))
  {
    return("mouse_cycle.rds")
  }
}

get_genome_annotation_names_mapping = function(spec, clean = F)
{
  if (spec %in% c("mouse","m"))
  {
    if(clean&&file.exists(annot_path("m")))
    {
      file.remove(annot_path("m"))
    }
    
    if(file.exists(annot_path("m")))
    {
      anno = readRDS(annot_path("m"))
    }
    else
    {
      if(testBioCConnection())
      {
        anno = getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name"), mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
        anno = list(gene_id = anno$ensembl_gene_id, gene_symbol = anno$mgi_symbol, chr = anno$chromosome_name)
        saveRDS(anno, file = annot_path("m"))
      }
      else
      {
        print("Internet connection required to load annotation for the first time")
        stop()
      }
    }
    
  }
  if (spec %in% c("human","h"))
  { 
    if(clean&&file.exists(annot_path("h")))
    {
      file.remove(annot_path("h"))
    }
    
    if(file.exists(annot_path("h")))
    {
      anno = readRDS(annot_path("h"))
    }
    else
    {
      if(testBioCConnection())
      {
        anno = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
        anno = list(ensembl_id = anno$ensembl_gene_id, gene_symbol = anno$hgnc_symbol, chr = anno$chromosome_name)
        saveRDS(anno, file = annot_path("h"))
      }
      else
      {
        print("Internet connection required to load annotation for the first time")
        stop()
      }
    }
  }
  return(anno)
}


detect_genome = function(gene_names)
{
  detect_genome_symbol_format = function(gene_names, annotation)
  {
    symbol = intersect(gene_names, annotation$gene_symbol)
    ensembl = intersect(gene_names, annotation$gene_id)
    
    if(length(symbol)  ==  0 && length(ensembl)  ==  0) return (NULL)
    if(length(symbol) > length(ensembl)) return (list(annot = "symbol", genes = symbol))
    else return (list(annot = "ensembl", genes = ensembl))
  }
  
  human = get_genome_annotation_names_mapping("h")
  mouse = get_genome_annotation_names_mapping("m")
  
  test_human = detect_genome_symbol_format(gene_names, human)
  test_mouse = detect_genome_symbol_format(gene_names, mouse)
  
  if(is.null(test_human)&&is.null(test_mouse))
  {
    print("Gene symbols should be ensembl or gencode. Either mouse or human")
    return(NULL)
  }
  
  if(!is.null(test_human))
  {
    if(is.null(test_mouse) || length(test_mouse$genes) < length(test_human$genes))
    {
      return(list(org = "human", annot = test_human$annot, genes = list(main = test_human$genes, collis = test_mouse$genes)))
    }
  }
  
  if(!is.null(test_mouse))
  {
    if(is.null(test_human) || length(test_mouse$genes) > length(test_human$genes))
    {
      return(list(org = "mouse", annot = test_mouse$annot, genes = list(main = test_mouse$genes, collis = test_human$genes)))
    }
  }
  
  # annot = NULL
  # if(!is.null(test_mouse)&&!is.null(test_mouse)&&test_mouse$annot  ==  test_human$annot)
  # {
  #   annot = test_mouse$annot
  # }
  # else
  # {
  #   if(!is.null(test_human)&&is.null(test_mouse))
  #   {
  #     annot = test_human$annot
  #   }
  # 
  #   if(!is.null(test_mouse)&&!is.null(test_human))
  #   {
  #     annot = test_mouse$annot
  #   }
  # }
  # 
  return(NULL)
}

get_annot_db = function(gene_names, org)
{
  db_to_use = NULL
  if(org %in% c("human", "h"))
  {
    db_to_use = "org.Hs.eg.db"
  }
  if(org %in% c("mouse", "m"))
  {
    db_to_use = "org.Mm.eg.db"
  }
  return(db_to_use)
}

gene_annot_convert_to = function(gene_names, new_encoding = "ENTREZID", org, annot)
{
  
  old_encoding = toupper(annot)
  
  org_db = get_annot_db(gene_names, org)
  
  new_encoding = toupper(new_encoding)
  
  if(new_encoding  ==  old_encoding)
  {
    return (list(old = gene_names, new = gene_names))
  }
  
  gene_annot_mapping = as.matrix(bitr(gene_names, fromType = old_encoding, toType = new_encoding, OrgDb = org_db))
  return(list(old = gene_annot_mapping[,1], new = gene_annot_mapping[,2]))
}

mtx_gene_annot_convert_to = function(mtx, new_encoding = "ENTREZID", org, annot)
{
  mtx = mtx[unique(rownames(mtx)),]
  annot_converted = gene_annot_convert_to(rownames(mtx), new_encoding = "ENSEMBL", org = org, annot = annot)
  names(annot_converted$new) = NULL
  
  mtx = mtx[match(annot_converted$old, rownames(mtx)),]
  rownames(mtx) = annot_converted$new
  
  mtx = mtx[unique(rownames(mtx)),]
  
  return(mtx)
}





# presort matrix by decreasing sum of row and columns
# can be usefull to select top-N expressed genes
sort_by_rowcol_sum = function(A)
{
  gene_sum = rowSums(A)
  cell_sum = c(colSums(A),0)
  
  B <-cbind(A, gene_sum)
  B <-rbind(B, cell_sum)
  
  genes_ordered = order(B[,"gene_sum"], decreasing = TRUE)
  cells_ordered = order(B["cell_sum",], decreasing = TRUE)
  
  B = B[genes_ordered, cells_ordered]
  
  sorted = B[-grep("cell_sum", rownames(B)), -grep("gene_sum", colnames(B))]
  return (sorted[rowSums(sorted)>0, colSums(sorted)>0])
}

rel2abs = function(mtx)
{
  rel2abs_base = function(mtx)
  {
    fd = as.matrix(rownames(mtx))
    rownames(fd) = fd
    colnames(fd)[1]<-"gene_short_name"
    
    pd = as.matrix(colnames(mtx))
    rownames(pd) = pd
    colnames(pd)[1]<-"cell_name"
    
    pd = new("AnnotatedDataFrame", data = as.data.frame(pd))
    fd = new("AnnotatedDataFrame", data = as.data.frame(fd))
    
    relative = newCellDataSet(as.matrix(mtx),
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = tobit(Lower = 0.1))
    
    rpc_matrix = sort_by_rowcol_sum(na.omit(t(na.omit(t(relative2abs(relative, t_estimate = estimate_t(exprs(relative)), method = "num_genes", cores = cores))))))
    
    return (rpc_matrix)
  }
  
  norm = tryCatch(rel2abs_base(mtx), error = function(e){"e"})
  if(typeof(norm) != "character")
  {
    return(norm) 
  }
  return(t(rel2abs_base(t(mtx))))
}


saver_impute <-function(expr, size_factors = NULL, genes = NULL, cores = parallel:::detectCores()-1)
{
  expr_new = sort_by_rowcol_sum(t(na.omit(t(na.omit(expr)))))
  if(is.null(genes))
  {
    genes = rownames(expr_new)
    if(nrow(expr_new)>3000)
    {
      genes = genes[1:3000]
    }
  }
  
  doParallel::stopImplicitCluster()
  doParallel::registerDoParallel(cores = cores)
  
  max_var = nrow(expr_new) 
  
  if(is.null(size_factors))
  {
    impute_sorted_in = saver(expr_new,
                             parallel = T, 
                             pred.genes = which(rownames(expr_new) %in% genes), 
                             nzero = max(10, ncol(expr_new)/20), 
                             dfmax = 5 #limit netglm model complexity to speed-up computations
                             #should not affect result too much, cuz every gene processed separately in SAVER
    )
  }
  else
  {
    impute_sorted_in = saver(expr_new,
                             parallel = T, 
                             size.factor = size_factors, 
                             pred.genes = which(rownames(expr_new) %in% genes), 
                             nzero = max(10, ncol(expr_new)/20),
                             dfmax = 5
    )
  }
  
  doParallel::stopImplicitCluster()
  
  rpc_matrix_new = t(na.omit(t(impute_sorted_in$estimate)))
  
  return (rpc_matrix_new)
}

matrix_sum = function(A, B)
{
  A[rownames(B), colnames(B)] = A[rownames(B), colnames(B)] + B[rownames(B), colnames(B)]
  return (A)
}

load_velocity_mtx = function(path)
{
  # Merge matrix of intron-exon matrix counts by gene name
  dropestOutput = readRDS(path)
  
  exo_mtx = as.matrix(dropestOutput$exon)
  intro_mtx = as.matrix(dropestOutput$intron)
  span_mtx = as.matrix(dropestOutput$spanning)
  
  detected_genes = unique(c(rownames(exo_mtx), rownames(intro_mtx), rownames(span_mtx)))
  cells_count = ncol(exo_mtx)
  expr_all = matrix(0, nrow = length(detected_genes), ncol = cells_count)
  
  rownames(expr_all) = detected_genes
  colnames(expr_all) = colnames(exo_mtx)
  
  expr_all = matrix_sum(expr_all, exo_mtx)
  expr_all = matrix_sum(expr_all, intro_mtx)
  expr_all = matrix_sum(expr_all, span_mtx)
  
  return(expr_all)
}

magic_impute = function (mtx, components = 20)
{
  magicimpute2matrix<-function(M)
  {
    B = data.frame(M)
    colnames(B) = colnames(M)
    rownames(B) = M[,"X"]
    B = t(as.matrix(B[,which(colnames(M)!= "X")]))
    storage.mode(B) = "numeric"
    return (B)
  }
  
  file_in = paste0(digest::digest(mtx), "_for_magic_exprs.csv")
  file_out = paste0(digest::digest(mtx), "_magic_out.csv")
  
  write.csv(mtx, file = file_in)
  
  tryCatch(system(paste0('bash -c "MAGIC.py -d ', file_in, ' -o ', file_out, '  --cell-axis \'columns\' -p ', min(components, nrow(mtx)),' csv"')))
  
  magic_imputed = tryCatch(na.omit(magicimpute2matrix(read.csv(file = file_out))), error = function(e){NULL})
  
  unlink(file_in)
  unlink(file_out)
  
  return (magic_imputed)
}


detect_DE_genes = function(mtx)
{
  Normalized_data = M3DropCleanData(mtx,
                                    labels = colnames(mtx),
                                    is.counts = TRUE)
  
  M3Drop_genes = M3DropGetExtremes(Normalized_data$data, fdr_threshold = 0.05)
  
  DE_genes = NULL
  for(method in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
  {
    DE_genes = c(DE_genes, rownames(M3DropDifferentialExpression(Normalized_data$data,
                                                                 mt_method = method, mt_threshold = 0.05, suppress.plot = T)))
  }
  
  DE_genes = unique(DE_genes, c(M3Drop_genes$left, M3Drop_genes$right))
  
  object = CreateSeuratObject(raw.data = mtx, min.cells = 3, min.genes = 200)
  object = NormalizeData(object = object, normalization.method = "LogNormalize", display.progress = F)
  object = FindVariableGenes(object = object, mean.function = ExpMean, dispersion.function = LogVMR)
  
  DE_genes = unique(c(object@var.genes, DE_genes))
  
  DE_genes = DE_genes[order(rowSums(mtx[DE_genes,]), decreasing = T)]
  
  return(DE_genes)
}



simple_preprocessing = function(sce, libsize = T, feature = T, mito = T, spikeIns = T)
{
  filter_mito = function (sce)
  {
    #MT QC
    detected_genome = detect_genome(rownames(sce))
    if(is.null(detected_genome))
    {
      return(NULL)
    }
    anno = get_genome_annotation_names_mapping(detected_genome$org)
    #Filter unknown genes
    is.mito =  rownames(sce) %in% anno$gene_symbol[which(anno$chr  ==  "MT")]
    
    is.mito_offline = (grepl("^mt-", rownames(sce)) | grepl("^MT-", rownames(sce)))
    
    if((exists("is.mito")&&length(which(is.mito)) == 0)||(!exists("is.mito")))
    {
      is.mito = is.mito_offline
    }
    return(rownames(sce)[which(is.mito)])
  }
  
  filter_spikeIns = function(sce)
  {
    is.spike = (grepl("^ERCC-", rownames(sce)) | grepl("^ercc-", rownames(sce)))
    isSpike(sce, "ercc") = is.spike
    return(rownames(sce)[which(is.spike)])
  }
  
  sce = calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(mt = filter_mito(sce), ercc = filter_spikeIns(sce)),
                           cell_controls = NULL, nmads = 3, pct_feature_controls_threshold = 80)
  
  final_drop = rep(F, ncol(sce))
  
  libsize.drop = rep(F, ncol(sce))
  if(libsize) {
    libsize.drop = isOutlier(sce$total_counts, nmads = 3, type = "low", log = T)
    if(!is.na(sum(libsize.drop)))
    {
      final_drop = libsize.drop|final_drop
    }
  }
  
  feature.drop = rep(F, ncol(sce))
  if(feature) {
    feature.drop = isOutlier(sce$total_features, nmads = 3, type = "low", log = T)
    if(!is.na(sum(feature.drop)))
    {
      final_drop = feature.drop|final_drop
    }
  }
  
  
  mito.drop = rep(F, ncol(sce))
  if(mito) {
    mito.drop = isOutlier(sce$total_counts_mt, nmads = 3, type = "high", log = T)
    if(!is.na(sum(mito.drop)))
    {
      final_drop = mito.drop|final_drop
    }
  }
  
  spike.drop = rep(F, ncol(sce))
  if(spikeIns) {
    spike.drop = isOutlier(sce$total_counts_ercc, nmads = 3, type = "both", log = T)
    if(!is.na(sum(spike.drop)))
    {
      final_drop = spike.drop|final_drop
    }
  }
  
  sce = sce[,!final_drop]
  print(data.frame(ByLibSize = sum(libsize.drop), ByFeature = sum(feature.drop),
                   ByMito = sum(mito.drop), BySpike = sum(spike.drop), Remaining = ncol(sce)))
  
  #Computing separate size factors for spike-in transcripts if exists
  sce1 = tryCatch(computeSpikeFactors(sce, type = "ercc", general.use = FALSE), error = function(e) {"e"})
  if(typeof(sce1) !=  "character")
  {
    sce = sce1
  }
  
  #Deconvolution method to normalize genes with lots of zero-counts
  sce = get_sum_factors(sce)
  #print(summary(sizeFactors(sce)))
  sce = normalise(sce)
  
  return(sce)
}

cell_cycle_markers = function(org)
{
  pairs = NULL
  if(org %in% c("human", "h"))
  {
    pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
  }
  
  if(org %in% c("mouse", "m"))
  {
    pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran"))
  }
  
  return(pairs)
}

cell_cycle_decomposed_normalization = function(sce)
{
  cyclone_mapping = function(sce, org, annot)
  {
    annot_map = get_genome_annotation_names_mapping(spec = org)
    if(annot  ==  "ensembl")
    {
      ensembl = annot_map$gene_id[match(rownames(sce), annot_map$gene_id)]
    }
    
    if(annot  ==  "symbol")
    {
      ensembl = annot_map$gene_id[match(rownames(sce), annot_map$gene_symbol)]
    }
    
    assignments = cyclone(sce, cell_cycle_markers(org), gene.names = ensembl)
    return (assignments)
  }
  
  detected_genome = detect_genome(rownames(sce))
  if(is.null(detected_genome$org) || is.null(detected_genome$annot))
  {
    print("Unable to deside which organism to work with, skipping cell_cycle_decomposition step")
  } else
  {
    assignments = cyclone_mapping(sce, detected_genome$org, detected_genome$annot)
    
    sce_g1 = sce[,which(assignments$normalized.scores$G1 >=  0.5)]
    sce_g2m = sce[,which(assignments$normalized.scores$G2M >=  0.5)]
    sce_s = sce[,which(assignments$normalized.scores$S >=  0.5)]
    
    return(list(G1 = exprs(simple_preprocessing(sce_g1)), G2M = exprs(simple_preprocessing(sce_g2m)), S = exprs(simple_preprocessing(sce_s))))
  }
  
  return("e")
}

saver_two_path_impute = function(sce)
{
  #detection of highly variable genes to reduce impute step required time
  #combine several methods to get best result after impute step
  #undetected DE genes should be passed unimputed to MAGIC impute
  DE_genes = tryCatch(detect_DE_genes(2^exprs(sce)-1), error = function(e){NULL})
  
  if(nrow(sce)>2000&&!is.null(DE_genes))
  {
    DE_genes_magic0 = DE_genes[seq(1, length(DE_genes),2)]
    length(DE_genes_magic0)
    imputed_genes0 = tryCatch(saver_impute(counts(sce), size_factors = sizeFactors(sce), genes = DE_genes_magic0), error = function(e){NULL})
    
    DE_genes_magic1 = DE_genes[seq(2, length(DE_genes),2)]
    length(DE_genes_magic1)
    imputed_genes1 = tryCatch(saver_impute(counts(sce), size_factors = sizeFactors(sce), genes = DE_genes_magic1), error = function(e){NULL})
    
    if(!is.null(imputed_genes0)&&!is.null(imputed_genes1))
    {
      return(list(x = exprs(simple_preprocessing(imputed_genes0)), y = exprs(simple_preprocessing(imputed_genes1))))
    }
  }else{
    imputed_genes = tryCatch(saver_impute(counts(sce), size_factors = sizeFactors(sce)), error = function(e){NULL})
    if(!is.null(imputed_genes))
    {
      return(exprs(simple_preprocessing(imputed_genes)))
    }else{
      print("Impute error, normalized count matrix returned")
      return(exprs(simple_preprocessing(sce)))
    }
  }
}


get_sum_factors = function(sce)
{
  min_l = max(2, ncol(sce)/100)
  max_l = min(ncol(sce)/10, 100)
  norm_levels = unique(round(seq(min_l, max_l, (max_l - min_l)/10)))
  
  sce = tryCatch(computeSumFactors(sce, sizes = norm_levels, positive = F),
                 warning = function(w)
                 {computeSumFactors(sce, sizes = norm_levels, positive = T)})
  return (sce)
}

sce_norm = function(mtx)
{
  sce_norm_base = function(mtx)
  {
    filter_simple = function(sce, numgenes_sh = 2, numcells_sh = 2, lowerDetectionLimit_gene = 1, lowerDetectionLimit_cell = 1)
    {
      #Simple quality control on the cells
      keep_feature_col = rowSums(counts(sce)>0) > 0
      keep_feature_row = colSums(counts(sce)>0) > 0
      sce = sce[keep_feature_col, keep_feature_row]
      
      numcells = nexprs(sce, lowerDetectionLimit = lowerDetectionLimit_cell, byrow = T)
      keep.gene = numcells >= numcells_sh
      
      numgenes = nexprs(sce, lowerDetectionLimit = lowerDetectionLimit_gene, byrow = F)
      keep.cell = numgenes >= numgenes_sh
      
      sce = sce[keep.gene, keep.cell]
      return(sce)
    }
    
    sce_temp = SingleCellExperiment(assays = list(counts = mtx), colData = get_types(colnames(mtx)))
    rowData(sce_temp)$feature_symbol = rownames(sce_temp)
    sce_temp = sce_temp[!duplicated(rowData(sce_temp)$feature_symbol), ]
    
    sce_temp = filter_simple(sce_temp)
    
    sce_temp = tryCatch(get_sum_factors(sce_temp))
    
    #print(summary(sizeFactors(sce_temp)))
    sce_temp<-normalise(sce_temp)
    return(sce_temp)
  }
  
  sce_temp = tryCatch(sce_norm_base(mtx), error = function(e){NULL})
  if(is.null(sce_temp))
  {
    sce_temp = sce_norm_base(rel2abs(mtx))
  }
  
  return (sce_temp)
}

full_preprocess = function(working_path = NULL, counts_mtx = NULL, cell_cycle = F, impute = F, clean = T)
{
  print("Loading data...")
  dir.create(working_path, showWarnings = FALSE)
  setwd(working_path)
  
  counts_input = NULL
  # Reuse final results if computed
  if(clean)
  {
    unlink("counts_final.rds")
    counts_input = counts_mtx
  }else{
    counts_input = loadCountsRds("counts_final.rds")
  }
  
  
  sce = sce_norm(counts_input)
  print(dim(sce))
  
  print("Simple preprocessing...")
  sce = simple_preprocessing(sce)
  print(dim(sce))
  
  if(!cell_cycle && !impute)
  {
    sce_saver = exprs(sce)
  }else{
    print("Cell cycle classification...")
    
    sce_cycled = "e"
    if(cell_cycle)
    {
      sce_cycled = tryCatch(cell_cycle_decomposed_normalization(sce), error = function(e) {"e"})
    }
    
    if(typeof(sce_cycled) !=  "character")
    {
      print(lapply(sce_cycled, dim))
      if(impute)
      {
        print("SI")
        sce_saver = lapply(sce_cycled, saver_two_path_impute)
      }else{
        print("S")
        sce_saver = lapply(sce_cycled, function(sce_set){return(exprs(simple_preprocessing(sce_norm(sce_set))))})
      }
    }
    else
    {
      if(impute)
      {
        print("I")
        sce_saver = saver_two_path_impute(sce)
      }else{
        print("N")
        sce_saver = exprs(sce)
      }
    }
  }
  print("Merge splited")
  sce_saver = 2^merge_splits(sce_saver)-1
  print("Done...")
  saveRDS(sce_saver, file = "counts_final.rds")
  return(sce_saver)
}




add_zero_genes = function (mtx, gene_list)
{
  new_extent = matrix(0L, nrow = length(gene_list), ncol = ncol(mtx))
  rownames(new_extent) = gene_list
  mtx = rbind(mtx, new_extent)
  return(mtx[order(rownames(mtx)),])
}

m_bind = function(A, B, lA = "", lB = "", lAe = "", lBe = "")
{ 
  genes_all = unique(c(rownames(A), rownames(B)))
  
  genes_zero_A = setdiff(genes_all, rownames(A))
  A = add_zero_genes(A, genes_zero_A)
  colnames(A) = sapply(1:ncol(A), function(x) {return(paste0(lA, unlist(strsplit(colnames(A)[x], ".", fixed = TRUE))[1], lAe, ".", x))})
  
  genes_zero_B = setdiff(genes_all, rownames(B))
  B = add_zero_genes(B, genes_zero_B)
  colnames(B) = sapply(1:ncol(B), function(x) {return(paste0(lB, unlist(strsplit(colnames(B)[x], ".", fixed = TRUE))[1], lBe, ".", x))})
  
  C = cbind(A[order(rownames(A)),], B[order(rownames(B)),])
  
  return (C)
}




get_types = function(arr)
{
  return(sapply(arr, function(cell) { return(unlist(strsplit(cell, ".", fixed = T))[1])}))
}


dub = function(dataset_ex, copy_times = 0, transposed = F)
{
  if(transposed)
  {
    dataset_ex = t(dataset_ex)
  }
  
  dataset_copy = dataset_ex
  temp_names = colnames(dataset_ex)
  if(copy_times == 0 || copy_times == 1)
  {
    return(dataset_ex)
  }
  
  for (i in 1:(copy_times-1))
  {
    colnames(dataset_copy) = sapply(temp_names, function(cell) { return(paste0(cell,".copy",i))})
    dataset_ex = cbind(dataset_ex, dataset_copy)
  }
  
  if(transposed)
  {
    return(t(dataset_ex))
  }else{
    return(dataset_ex)
  }
}

get_cell_names_by_ident = function(object, clusters)
{
  ident = get_ident(object)
  
  return(names(ident)[which(ident%in%clusters)])
}

get_ident = function(object, types = NULL)
{
  ident = object@ident[which(!grepl(".copy", names(object@ident)))]
  
  if(!is.null(types))
  {
    return(ident[which(get_types(names(ident))%in%types)])
  }
  return(ident)
} 

get_orig_ident = function(object, types = NULL)
{
  ident = object@meta.data$orig.ident[which(!grepl(".copy", names(object@ident)))]
  
  if(!is.null(types))
  {
    return(ident[which(get_types(names(ident))%in%types)])
  }
  return(ident)
} 






merge_saver_split = function(saver_splited, lA = ".x", lB = ".y", lC = ".z")
{
  A = saver_splited$x
  B = saver_splited$y
  
  genes_all = unique(c(rownames(A), rownames(B), rownames(C)))
  
  genes_zero_A = setdiff(genes_all, rownames(A))
  A = add_zero_genes(A, genes_zero_A)
  colnames(A) = sapply(1:ncol(A), function(x) {return(paste0(unlist(strsplit(colnames(A)[x], ".", fixed = TRUE))[1], lA, ".", x))})
  
  genes_zero_B = setdiff(genes_all, rownames(B))
  B = add_zero_genes(B, genes_zero_B)
  colnames(B) = sapply(1:ncol(B), function(x) {return(paste0(unlist(strsplit(colnames(B)[x], ".", fixed = TRUE))[1], lB, ".", x))})
  
  E = cbind(A[order(rownames(A)),], B[order(rownames(B)),])
  
  return (E)
}

merge_cycle = function(cycle_splited, lA = ".g1", lB = ".g2m", lC = ".s")
{
  if("x" %in% names(cycle_splited$G1))
  {
    A = merge_saver_split(cycle_splited$G1)
    B = merge_saver_split(cycle_splited$G2M)
    C = merge_saver_split(cycle_splited$S)
  }else{
    A = cycle_splited$G1
    B = cycle_splited$G2M
    C = cycle_splited$S
  }
  genes_all = unique(c(rownames(A), rownames(B), rownames(C)))
  
  genes_zero_A = setdiff(genes_all, rownames(A))
  A = add_zero_genes(A, genes_zero_A)
  colnames(A) = sapply(1:ncol(A), function(x) {return(paste0(unlist(strsplit(colnames(A)[x], ".", fixed = TRUE))[1], lA, ".", x))})
  
  genes_zero_B = setdiff(genes_all, rownames(B))
  B = add_zero_genes(B, genes_zero_B)
  colnames(B) = sapply(1:ncol(B), function(x) {return(paste0(unlist(strsplit(colnames(B)[x], ".", fixed = TRUE))[1], lB, ".", x))})
  
  genes_zero_C = setdiff(genes_all, rownames(C))
  C = add_zero_genes(C, genes_zero_C)
  colnames(C) = sapply(1:ncol(C), function(x) {return(paste0(unlist(strsplit(colnames(C)[x], ".", fixed = TRUE))[1], lC, ".", x))})
  
  E = cbind(A[order(rownames(A)),], B[order(rownames(B)),], C[order(rownames(C)),])
  return(E)
}

merge_splits = function(splits)
{
  if("G1" %in% names(splits) && !is.null(splits$G1))
  {
    return(as(merge_cycle(splits), "dgCMatrix"))
  }
  if("x" %in% names(splits) &&!is.null(splits$x))
  {
    return(as(merge_saver_split(splits), "dgCMatrix"))
  }
  return(na.omit(splits))
}







try_assign_cell_cycle = function(object)
{
  genome = detect_genome(rownames(object@raw.data))
  if(is.null(genome))
  {
    return(NULL)
  }
  cc = cell_cycle_markers(genome$org)
  G2M = gene_annot_convert_to(unique(cc$G2M$first), new_encoding = genome$annot, org = genome$org, annot = "ENSEMBL")$new
  S = gene_annot_convert_to(unique(cc$S$first), new_encoding = genome$annot, org = genome$org, annot = "ENSEMBL")$new
  object1 = tryCatch(CellCycleScoring(object, G2M, S, set.ident = F), error = function(e){NULL})
  return(object1)
}

seur_norm_regress = function(object, regress, regress_cc, nbin, cutoff)
{
  object@meta.data[, "orig.ident"] = get_types(colnames(object@raw.data))
  object = NormalizeData(object = object, normalization.method = NULL, display.progress = F)
  vars_to_regress = NULL
  if(regress)
  {
    vars_to_regress = c(vars_to_regress, "nUMI")
  }
  
  if(!is.null(regress_cc))
  {
    object1 = try_assign_cell_cycle(object)
    if(!is.null(object1))
    {
      object = object1
      rm(object1)
      if(regress_cc == "cc_diff")
      {
        object@meta.data$CC.Difference = object@meta.data$S.Score - object@meta.data$G2M.Score
        vars_to_regress = c(vars_to_regress, "CC.Difference")
      }
      else{
        vars_to_regress = c(vars_to_regress, "S.Score", "G2M.Score")
      }
    }
  }
  
  if(!is.null(vars_to_regress))
  {
    object = ScaleData(object = object, vars.to.regress = vars_to_regress, display.progress = F, check.for.norm = F)
  }else{
    object = ScaleData(object = object, display.progress = F, check.for.norm = F)
  }
  
  object = FindVariableGenes(object, display.progress = F, do.plot = F, x.low.cutoff = -Inf, x.high.cutoff = Inf, y.cutoff = cutoff, num.bin = nbin)
  
  return(object)
}

compute_pca = function(object)
{
  pc_seq = c(80, 70, 60, 50, 40, 30, 20, 15, 12, 10, 7, 6, 5, 4, 3)
  
  
  for(pc_use in pc_seq[pc_seq<ncol(object@raw.data)])
  {
    object1 = tryCatch(RunPCA(object = object, pc.genes = object@var.genes#rownames(object@raw.data)
                              , print.results = FALSE, pcs.print = NULL, do.print = F, pcs.compute = pc_use, rev.pca = F, weight.by.var = T),
                       error = function(e) {"e"}, 
                       warning = function(w) {"w"}
    )
    
    if(typeof(object1) !=  "character")
    {
      #print(paste0("Pc_use: ", pc_use))
      # print(paste0("Determine statistical significance of PCA scores..."))
      # object1 <- JackStraw(object1, num.pc = ncol(object1@dr$pca@cell.embeddings), num.replicate = 100, prop.freq = 0.1)
      return (object1)
    }
  }
  return (NULL)
}

select_good_PCs = function(object)
{
  
  calc_pval = function(object, p_val = 1e-5)
  {
    pAll <- object@dr$pca@jackstraw@emperical.p.value
    pAll <- pAll[, 1:ncol(object@dr$pca@cell.embeddings), drop = FALSE]
    pAll <- as.data.frame(pAll)
    
    pcs_all = rep(1, ncol(object@dr$pca@cell.embeddings))
    names(pcs_all) = 1:ncol(object@dr$pca@cell.embeddings)
    
    for (i in 1:ncol(object@dr$pca@cell.embeddings))
    {
      pc.score <- suppressWarnings(prop.test(
        x = c(
          length(which(pAll[, i] <= p_val)),
          floor(nrow(pAll) * p_val)
        ),
        n = c(nrow(pAll), nrow(pAll))
      )$p.val)
      if (length(which(pAll[, i] <= p_val)) == 0) {
        pc.score <- 1
      }
      pcs_all[i] = pc.score
    }
    
    pcs_all = pcs_all[order(pcs_all)][which(pcs_all <= p_val)]
    
    return(pcs_all)
  }
  
  p_val = 1e-2
  pc_good = calc_pval(object, p_val = p_val)
  # while(length(pc_good)>=5)
  # {
  #   if(p_val<=1e-8)
  #   {
  #     break
  #   }
  #   p_val = p_val/10
  #   pc_good = calc_pval(object, p_val = p_val)
  #   print(pc_good)
  # }
  # 
  # while(length(pc_good)<=4)
  # {
  #   if(p_val == 1)
  #   {
  #     break
  #   }
  #   p_val = p_val*10
  #   pc_good = calc_pval(object, p_val = p_val)
  print(pc_good)
  # }
  
  cat("PC_pval: ", p_val, "\n")
  
  # if(p_val == 1)
  # {
  #   return(NULL)
  # }else{
  return(sort(as.numeric(names(pc_good[1:5]))))
  # }
}

compute_tsne = function(object, dim.embed = 2, reduction.use = "pca")
{
  init_plex = max(10, ncol(object@raw.data)/4)
  end_plex = 2
  step = (init_plex-end_plex)/5
  
  
  if(is.null(reduction.use))
  {
    genes.use = select_good_genes(object)
  }else{
    dim.use = 1:ncol(object@dr[[reduction.use]]@cell.embeddings)
  }
  
  # if(reduction.use == "pca")
  # {
  #   dim.use = select_good_PCs(object)
  # }
  
  
  for(perplex in unique(round(seq(init_plex, end_plex, -step))))
  {
    
    if(!is.null(reduction.use))
    {
      object1 = tryCatch(RunTSNE(object = object,
                                 dims.use = dim.use, 
                                 reduction.use = reduction.use,
                                 dim.embed = dim.embed, check_duplicates = FALSE, perplexity = perplex),error = function(e) {NULL})
    }else{
      
      object1 = tryCatch(RunTSNE(object = object, 
                                 genes.use = genes.use,
                                 dim.embed = dim.embed, check_duplicates = FALSE, perplexity = perplex),error = function(e) {NULL})
    }
    
    
    if(!is.null(object1))
    {
      #print(paste0("perplex ", perplex))
      return(object1)
    }
  }
  return(NULL)
}


compute_clustering = function(object, method = "min")
{
  #print("Compute tsne...")
  
  object = compute_tsne(object, dim.embed = 2, reduction.use = "pca")
  #good_PC = select_good_PCs(object)
  if(is.null(object)){
    return(NULL)
  }
  k.scale = 2
  
  #print("Clustering...")
  test_clustering = function(object, k.param = 10, k.scale = 2, target_clusters = 2)
  {
    test_cl = function(clusters)
    {
      if(length(clusters)>target_clusters)
      {
        return(2)
      }else if(length(clusters) < target_clusters)
      {
        return(1)
      }else if(length(clusters) == target_clusters)
      {
        if(min(clusters)>k.scale)
        {
          return(0)
        }else{
          return(1)
        }
      }
    }
    
    s_out <- tryCatch(capture.output(
      {object = FindClusters(object = object,
                             dims.use = 1:ncol(object@dr$tsne@cell.embeddings),
                             reduction.type = "tsne", 
                             #resolution = 0, 
                             k.param = k.param, k.scale = k.scale, algorithm = 3,
                             modularity.fxn = 1, print.output = F, save.SNN = F, force.recalc = F, 
                             #prune.SNN = 0, 
                             n.start = 10)}
    ), error = function(e){NULL})
    
    s_out <- unlist(strsplit(s_out[which(grepl("singletons identified", s_out))], c(" |\""), fixed = F))
    s_out <- s_out[which(grepl("^\\d$", s_out))][1]
    
    if(!is.null(s_out))
    {
      return(list(obj = object, test = 2))
    }else{
      return(list(obj = object, test = test_cl(as.vector(table(get_ident(object))))))
    }
  }
  
  good_rez = NULL
  direction = 1
  max_num = 1
  from = k.scale+1
  to = from+100
  test_rez = NULL
  
  for(depth in max_num:0)
  {
    if(from==to)
    {
      break
    }
    
    for(rez in seq(from, to, direction*(10^depth)))
    {
      object = test_clustering(object, k.param = rez)
      test_rez = object$test
      object = object$obj
      
      # print(rez)
      # print(test_rez)
      
      if(method == "min")
      {
        if(direction > 0)
        {
          if (test_rez == 2) {
            next
          }else if (test_rez == 0) {
            good_rez = rez
            to = from
            from = rez
            break
          }else if(test_rez == 1){
            good_rez = rez
            to = from
            from = rez
            break
          }
        }else{
          if (test_rez == 2) {
            to = from
            from = rez
            break
          }else if (test_rez == 0) {
            good_rez = rez
          }
        }
      }else if(method == "max")
      {
        if(direction > 0)
        {
          if (test_rez == 2) {
            next
          }else if (test_rez == 0) {
            good_rez = rez
            from = rez
            #maintain the direction
            direction = -direction
            break
          }else if(test_rez == 1){
            good_rez = rez
            to = from
            from = rez
            break
          }
        }else{
          if (test_rez == 2) {
            to = from
            from = rez
            break
          }else if (test_rez == 0) {
            good_rez = rez
            break
          }
        }
      }
      
    }
    direction = -direction
  }
  
  if(is.null(good_rez))
  {
    return(NULL)
  }
  
  object = test_clustering(object, k.param = good_rez)
  
  return(object$obj)
}

compute_clustering_min = function(object, strict_binary = T)
{
  #print("Find optimal PCs...")
  #good_PC = select_good_PCs(object)
  #print(good_PC)
  
  print("Compute tsne...")
  
  object = compute_tsne(object, dim.embed = 2)
  
  k.param_factor = 2
  k.param = max(3, ncol(object@raw.data)/k.param_factor)
  k.scale = 2
  gran_shresh = k.scale
  
  print("Building SNN...")
  object = FindClusters(object = object,
                        dims.use = 1:ncol(object@dr$tsne@cell.embeddings)
                        , reduction.type = "tsne", 
                        resolution = 0, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T, prune.SNN = 0, n.start = 10)
  
  print("Clustering...")
  test_clustering = function(clusters)
  {
    if(length(clusters)<2)
    {
      return(2)
    }else if(length(clusters)>= 2)
    {
      if(length(clusters) == 2)
      {
        if(min(clusters)>gran_shresh)
        {
          return(0)
        }else{
          return(2)
        }
      }else if(length(clusters) > 2)
      {
        if(length(which(clusters>gran_shresh)>1))
        {
          return(1)
        }else{
          return(NULL) 
        }
      }
    }else{
      return(-1)
    }
  }
  
  good_rez = NULL
  direction = 1
  max_num = 10
  from = 0
  to = 30
  test_rez = NULL
  
  last_clusters = NULL
  clusters = NULL
  last_count = 1
  
  for(depth in 0:max_num)
  {
    if(from==to)
    {
      break
    }
    
    for(rez in seq(from, to, direction*(10^-depth)))
    {
      if(!strict_binary)
      {
        if(!is.null(clusters) && !is.null(last_clusters) && all(as.vector(last_clusters)==as.vector(clusters)))
        {
          last_count = last_count+1
        }else{
          last_count = 1
          last_clusters = clusters
        }
      }
      
      object = FindClusters(object = object, resolution = rez, algorithm = 3, modularity.fxn = 1, print.output = F, n.start = 10, reuse.SNN = T)
      
      clusters = as.vector(table(get_ident(object)))
      test_rez = test_clustering(clusters)
      
      #print(test_rez)
      if(direction > 0)
      {
        if(is.null(test_rez)){
          return(NULL)
        }else if (test_rez == 2) {
          next
        }else if (test_rez == 0) {
          good_rez = rez
          break
        }else if(test_rez == 1){
          good_rez = rez
          to = from
          from = rez
          break
        }else if(test_rez == -1){
          to = from
          from = rez
          break
        }
      }else{
        if(is.null(test_rez)){
          next
        }else if (test_rez == 2) {
          to = from
          from = rez
          break
        }else if (test_rez == 0) {
          good_rez = rez
          break
        }else if(test_rez == 1){
          good_rez = rez
          next
        }else if(test_rez == -1){
          next
        }
      }
    }
    
    if(test_rez == 0 || (!strict_binary&&test_rez == 1))
    {
      break
    }
    direction = -direction
  }
  
  object = FindClusters(object = object, resolution = good_rez, algorithm = 3, modularity.fxn = 1, print.output = F, n.start = 10, reuse.SNN = T)
  
  return(object)
}


seurat_analyse_mtx = function(object, regress, regress_cc, do.magic, do.cluster = T)
{
  raw = object
  object = exprs(sce_norm(raw))
  
  if(do.magic)
  {
    magic = magic_impute(as.matrix(object))
    if(!is.null(magic))
    {
      print("magic")
      raw = as(raw[,colnames(magic)], "dgCMatrix")
      object = magic
    }
  }
  
  object = CreateSeuratObject(raw.data = as(object, "dgCMatrix"))
  
  #
  # for(nbin in seq(10, 50, 10))
  # {
  #   for(cutoff in c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1))
  #   {
  #     object1 = tryCatch(seurat_analyse(object, regress = regress, regress_cc = regress_cc, do.cluster = do.cluster, raw_before_magic = raw, nbin = round(ncol(object@raw.data)/nbin), cutoff = cutoff),error = function(e){NULL})
  #     if(!is.null(object1)){
  #     if(all(rowMins(table(get_orig_ident(object1), get_ident(object1)))!=0))
  #       {
  #       cat("\n mult: ", nbin," nbin: ", round(ncol(object@raw.data)/nbin), " cutoff: ", cutoff," F", "\n", sep = "\t")
  #     }else if(all(rowMins(table(get_orig_ident(object1), get_ident(object1)))==0)){
  #       cat("\n mult: ", nbin," nbin: ", round(ncol(object1@raw.data)/nbin), " cutoff: ", cutoff," T", "\n", sep = "\t")
  #     }else if(any(rowMins(table(get_orig_ident(object1), get_ident(object1)))==0)){
  #       cat("\n mult: ", nbin," nbin: ", round(ncol(object1@raw.data)/nbin), " cutoff: ", cutoff," TP", "\n", sep = "\t")
  #     }
  #     }
  #   }
  # }
  
  return(seurat_analyse(object, regress = regress, regress_cc = regress_cc, do.cluster = do.cluster, raw_before_magic = raw))
}

seurat_analyse = function(object, regress, regress_cc, do.cluster = T, raw_before_magic = NULL, nbin = 20, cutoff = 0.2)
{
  object = seur_norm_regress(object, regress = regress, regress_cc = regress_cc, nbin = nbin, cutoff = cutoff)
  
  object1 = compute_pca(object)
  if(!is.null(object1))
  {
    object = object1
    rm(object1)
  }else{
    return("e")
  }
  
  if(do.cluster)
  {
    object1 = compute_clustering(object, method = "min")
    if(!is.null(object1))
    {
      object = object1
      rm(object1)
    }else{
      return("e")
    }
  }
  if(!is.null(raw_before_magic))
  {
    object@raw.data <- raw_before_magic
  }
  
  return (object)
}

get_copy_count = function(object)
{
  if("raw.data" %in% slotNames(object))
  {
    object = object@raw.data
  }
  
  first_split = sapply(colnames(object), function(cell) { return(unlist(strsplit(cell, ".copy", fixed = TRUE))[2])})
  mult = max(as.numeric(first_split[which(!is.na(first_split))]))+1
  if(is.infinite(mult))
  {
    mult = 1
  }
  return(mult)
}

rename_clusters = function(object, parent.ident = NULL)
{
  if(is.null(parent.ident))
  {
    parent.ident = "s"
  }
  object_cl = GetClusters(object)
  object_cl$cluster = match(object_cl$cluster, unique(object_cl$cluster))
  object_cl$cluster = lapply(object_cl$cluster, function(new_ident){return(paste0(parent.ident, new_ident))})
  object = SetClusters(object, object_cl)
  return(object)
}

pick_last_char = function(x)
{
  return(substr(x, nchar(x), nchar(x)))
}

seur_get_raw_and_params = function(object, ident = NULL, regress, regress_cc, do.magic, search_markers)
{
  if(!is.null(ident) && ident  == "SeuratProject"){
    ident = NULL
  }
  #packing params for parallel run
  if(is.null(ident))
  {
    raw_data = object@raw.data
  }else{
    raw_data = object@raw.data[, get_cell_names_by_ident(object, ident)]
  }
  return(list(object = raw_data, ident = ident, regress = regress, regress_cc = regress_cc, do.magic = do.magic, search_markers = search_markers))
}

finalize_no_split = function(object, ident)
{
  object = CreateSeuratObject(raw.data = as(object, "dgCMatrix"))
  
  cl_labels = GetClusters(object)
  if(pick_last_char(toString(ident)) !=  "f")
  {
    cl_labels$cluster = rep(paste0(ident, "f"), length(cl_labels$cluster))
  }else{
    cl_labels$cluster = rep(ident, length(cl_labels$cluster))
  }
  object = SetClusters(object, cl_labels)
  
  return(object)
}

seur_find_split = function(object)
{
  ident = object$ident
  regress = object$regress
  regress_cc = object$regress_cc
  do.magic = object$do.magic
  raw_data = object$object
  search_markers = object$search_markers
  
  #curr_mult = get_copy_count(raw_data)
  if(is.null(ident) || pick_last_char(toString(ident)) !=  "f")
  {
    object = tryCatch(seurat_analyse_mtx(raw_data, regress = regress, regress_cc = regress_cc, do.magic = do.magic, do.cluster = T), error = function(e) {print(e)
      return("e")})
    if(typeof(object) !=  "character")
    {
      if(length(table(get_ident(object))) == ncol(raw_data))
      {
        break
      }
      
      if(length(table(get_ident(object)))>1)
      {
        object = rename_clusters(object, ident)
        cl_markers = list()
        if(search_markers)
        {
          cl_markers = tryCatch(FindAllMarkers(object, test.use = "MAST", logfc.threshold = 0.05, return.thresh = 0.01, only.pos = F, min.pct = 0.7), error = function(e){NULL})
          if(length(cl_markers)==0)
          {
            object = finalize_no_split(raw_data, ident)
            return(list(obj = object, mark = list()))
          }
        }
        
        return(list(obj = object, mark = cl_markers))
      }
    }
  }
  
  object = finalize_no_split(raw_data, ident)
  
  return(list(obj = object, mark = list()))
}


seur_stat = function(object)
{
  table(get_ident(object), get_orig_ident(object))
  
  plot_seur(object)
  
  plot_seur(object, method =  "pca")
}

split_ident = function(object, ident = NULL, regress = F, regress_cc = NULL, do.magic = F, search_markers = T)
{
  param_list_for_ident = seur_get_raw_and_params(object, ident = ident, regress = regress, regress_cc = regress_cc, do.magic = do.magic, search_markers = search_markers)
  return(seur_find_split(param_list_for_ident))
}

split_all_ident = function(object, regress = F, regress_cc = NULL, do.magic = F, search_markers = T)
{
  parallel_list = lapply(unique(get_ident(object)), function(ident){ seur_get_raw_and_params(object, ident = ident, regress = regress, regress_cc =  regress_cc, do.magic = do.magic, search_markers = search_markers)})
  
  #clean threads if the last run was terminated incorrectly
  tryCatch(doParallel::stopImplicitCluster(),error = function(e){NULL})
  
  doParallel::registerDoParallel(cores = parallel::detectCores()-1)
  
  parallel_list = foreach(object = parallel_list, .options.snow = list(preschedule = F)) %dopar% {
    source('Counts2Exprs.R', local = TRUE)
    object = seur_find_split(object)
    return(object)
  }
  
  doParallel::stopImplicitCluster()
  
  return(parallel_list)
}

seur_merge = function(object)
{
  if(!is.null(object$obj)){
    new_idents = GetClusters(object$obj)
    merged_object = object$obj
    mark = object$mark
  }else
  {
    cl_idents = lapply(1:length(object), function(i) {return(GetClusters(object[[i]]$obj))})
    new_idents = cl_idents[[1]]
    merged_object = object[[1]]$obj
    mark = object[[1]]$mark
    
    object[[1]] = "NULL"
    
    for(i in 2:length(cl_idents))
    {
      new_idents = rbind(new_idents, cl_idents[[i]])
      merged_object = MergeSeurat(merged_object, object[[i]]$obj, do.normalize = F)
      mark = rbind(mark, object[[i]]$mark)
      
      object[[i]] = "NULL"
    }
  }
  merged_object = SetClusters(merged_object, new_idents)
  merged_object = seurat_analyse(merged_object, regress = F, regress_cc = NULL, do.cluster = F)
  
  return (list(obj = merged_object, mark = mark))
}





plot_seur<- function(object, method = "tsne")
{
  if(method  ==  "tsne")
  {
    if(is.null(object@dr$tsne))
    {
      object = compute_tsne(object)
    }
    p1 = TSNEPlot(object = object, group.by = "orig.ident", do.return = TRUE, pt.size = 1.5)
    p2 = TSNEPlot(object = object, do.return = TRUE, pt.size = 1.5)
    plot_grid(p1, p2)
  }else
    if(method  ==  "pca")
    {
      p1 = PCAPlot(object = object, group.by = "orig.ident", dim.1 = 1, dim.2 = 2, do.return = TRUE, pt.size = 1.5)
      p2 = PCAPlot(object = object, do.return = TRUE, dim.1 = 1, dim.2 = 2, pt.size = 1.5)
      plot_grid(p1, p2)
    }
}

sort_diag = function(object)
{
  return(object[, order(sapply(1:ncol(object), function (x) {return(which.max(object[,x]))}))])
}


cluster_recursive = function(object, regress = F, regress_cc = NULL, do.magic = F, lvls = NULL, search_markers = T)
{
  #non-parallel on fist split
  object = split_ident(CreateSeuratObject(raw.data = object), regress = regress, regress_cc = regress_cc, do.magic = do.magic, search_markers = search_markers)
  print(sort_diag(table(get_orig_ident(object$obj), get_ident(object$obj))))
  final_markers = object$mark
  
  while(any(sapply(sapply(get_ident(object$obj), toString), pick_last_char)!="f") || (!is.null(lvls)&&lvls>1))
  {
    if(!is.null(lvls))
    {
      lvls = lvls-1
    }
    
    object = split_all_ident(object$obj, regress = regress, regress_cc = regress_cc, do.magic = do.magic, search_markers = search_markers)
    object = seur_merge(object)
    
    final_markers = rbind(final_markers, object$mark)
    
    print(sort_diag(table(get_orig_ident(object$obj), get_ident(object$obj))))
  }
  
  return(list(obj = object$obj, mark = final_markers))
}



plot_seur_3d = function(object, method = "tsne", radius = 0.3, lvls = NULL, old_ident = NULL)
{
  library(plyr)
  library(scatterplot3d)
  library(rgl)
  
  if(method  ==  "tsne")
  {
    object = compute_tsne(object, dim.embed = 3)
    rl = cbind(
      object@dr$tsne@cell.embeddings[,1],
      object@dr$tsne@cell.embeddings[,2],
      object@dr$tsne@cell.embeddings[,3]
    )
  }
  
  if(method  ==  "pca")
  {
    rl = cbind(
      object@dr$pca@cell.embeddings[,1],
      object@dr$pca@cell.embeddings[,2],
      object@dr$pca@cell.embeddings[,3]
    ) 
    
  }
  
  if(is.null(lvls))
  {
    curr_labels = get_ident(object)
  }else{
    curr_labels = substr(sapply(get_ident(object),toString), start = 1, stop = lvls)
  }
  
  cells_lab = unique(curr_labels)
  labels_uniq = data.frame(names = cells_lab, colors =  rainbow(length(cells_lab)))
  
  mfrow3d(1, 2)
  next3d()
  plot3d(x = rl[,1], y = rl[,2], z = rl[,3], col = mapvalues(curr_labels, from = array(labels_uniq$names), to = array(labels_uniq$colors)), type = "s", radius = radius)
  legend3d("topright", col = array(labels_uniq$colors), legend = labels_uniq$names, cex = 2, pch = 0.1)
  next3d()
  
  if(is.null(old_ident))
  {
    curr_labels = get_orig_ident(object)
  }else{
    curr_labels = old_ident
  }
  
  cells_lab = unique(curr_labels)
  labels_uniq = data.frame(names = cells_lab, colors =  rainbow(length(cells_lab)))
  
  plot3d(x = rl[,1], y = rl[,2], z = rl[,3], col = mapvalues(curr_labels, from = array(labels_uniq$names), to = array(labels_uniq$colors)), type = "s", radius = radius)
  legend3d("topright", col = array(labels_uniq$colors), legend = labels_uniq$names, cex = 2, pch = 0.1)  
  
}



#GO
# plot_venn = function(cl_markers)
# {
#   VENN.LIST = cl_markers$cluster.markers.list[1:min(5, length(cl_markers$cluster.markers))]
#   VENN.LIST = lapply(VENN.LIST, sort)
# 
#   names(VENN.LIST) = (1:min(5, length(cl_markers$cluster.markers)))-1
#   grid.newpage()
#   venn.plot = venn.diagram(VENN.LIST, NULL, fill = rainbow(length(names(VENN.LIST))), 
#                              alpha = rep(0.5,length(names(VENN.LIST))), cex = 2, cat.fontface = 4, 
#                              category.names = names(VENN.LIST), main = "DE Genes Inters")
#  
#   grid.draw(venn.plot)
# }

# find_markers = function(object, p_val = 0.05, test.use = "MAST", logfc.threshold = 0.5, only.pos = T, min.pct = 0.8, min.diff.pct = -Inf)
# {
#   cluster.markers = NULL
#   cluster_numbers = NULL
#   
#   for (ident in sort(unique(get_ident(object))))
#   {
#     fm = tryCatch(FindMarkers(object = object, ident.1 = ident,  test.use = test.use, logfc.threshold = logfc.threshold, only.pos = only.pos, min.pct = min.pct, min.diff.pct = min.diff.pct), error = function(e){NULL})
#     if(!is.null(fm))
#     {
#       cluster.markers = c(cluster.markers, list(ident = fm))
#       cluster_numbers = c(cluster_numbers, ident)
#     }
#   }
#   
#   #filter by p_val
#   for (i in 1:length(cluster.markers))
#   {
#     if(length(cluster.markers[i]$ident) == 0){ next }
#     
#     cluster.markers[i]$ident = cluster.markers[i]$ident[which(cluster.markers[i]$ident$p_val_adj<= p_val),]
#   }
#   
#   
#   cluster.markers.list = NULL
#   for (i in 1:length(cluster.markers))
#   {
#     if(length(cluster.markers[i]$ident) == 0){
#       cluster.markers.list = c(cluster.markers.list, list(item = list()))
#       next }
# 
#     genome<-detect_genome(rownames(cluster.markers[i]$ident))
#     if(!is.null(genome))
#     {
#       entrez_converted = gene_annot_convert_to(rownames(cluster.markers[i]$ident), org = genome$org, annot = genome$annot)
#       names(entrez_converted$new) = NULL
#   
#       cluster.markers[i]$ident = cluster.markers[i]$ident[match(entrez_converted$old, rownames(cluster.markers[i]$ident)),]
#       cluster.markers.list = c(cluster.markers.list, list(item = entrez_converted$new))
#     
#     }
#     cluster.markers[i]$ident = cluster.markers[i]$ident[order(cluster.markers[i]$ident$avg_logFC, decreasing = T),]
#   }
#   cl_markers = list(cluster.markers = cluster.markers, cluster.markers.list = cluster.markers.list)
#   
#   return(cl_markers)
# }
# 
# GO_enr = function(cl_markers, p_val = 0.05, q_val = 0.05, new_window = T, org = "m")
# {
# 
#   cluster.markers = cl_markers$cluster.markers
#   cluster.markers.list = cl_markers$cluster.markers.list
# 
#   names(cluster.markers.list) = rep("item", length(cluster.markers.list))
#   for (i in 1:length(cluster.markers.list))
#   {
#     cl_mark = cluster.markers.list[i]$item
# 
#     x = enrichGO(gene = cl_mark, pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T, OrgDb = get_annot_db(rownames(cluster.markers$ident), org = org))
#     if(nrow(x@result)>0)
#     {
#       x@result$Description[order(x@result$GeneRatio,decreasing = T)]
#       enrichMap(x, layout = igraph::layout_components, vertex.label.cex = 1, fixed = !new_window)
#       cnetplot(x, categorySize = "GeneRatio", fixed = !new_window)
#     }
#   }
#   names(cluster.markers.list) =  sapply(1:length(cluster.markers.list), function(cl){return(paste0("X",cl-1))})
# 
#   res = tryCatch(compareCluster(cluster.markers.list, fun = "enrichGO", OrgDb = get_annot_db(rownames(cluster.markers$ident))), error = function(e){"e"})
#     if(typeof(res) !=  "character")
#     {
#       plot(res, showCategory = 20)
#     }
# }


