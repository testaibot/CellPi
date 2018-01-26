


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

magic_impute = function (mtx, components = 100)
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
  object = NormalizeData(object = object, normalization.method = "LogNormalize")
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
  
  if(any(mtx%%1!=0))
  {
    mtx = rel2abs(log2(mtx+1))
  }
  
  sce_temp = SingleCellExperiment(assays = list(counts = mtx), colData = get_types(colnames(mtx)))
  rowData(sce_temp)$feature_symbol = rownames(sce_temp)
  sce_temp = sce_temp[!duplicated(rowData(sce_temp)$feature_symbol), ]
  
  sce_temp = filter_simple(sce_temp)
  
  sce_temp = get_sum_factors(sce_temp)
  
  #print(summary(sizeFactors(sce_temp)))
  sce_temp<-normalize(sce_temp)
  return (sce_temp)
}

full_preprocess = function(working_path = NULL, counts_mtx = NULL, cell_cycle = T, impute = F, markers = NULL, clean = FALSE)
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

seur_norm_regress = function(object, regress, regress_cc)
{
  object@meta.data[, "orig.ident"] = get_types(colnames(object@raw.data))
  object = NormalizeData(object = object, normalization.method = "LogNormalize")
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
      object<-object1
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
    object = ScaleData(object = object, vars.to.regress = vars_to_regress, display.progress = T)
  }else{
    object = ScaleData(object = object, display.progress = T)
  }
  
  
  return(object)
}

compute_pca = function(object)
{
  init_pc = min(max(nrow(object@raw.data)/120, 30), 100)
  end_pc = 2
  step = (init_pc-end_pc)/10
  
  for(pc_use in unique(round(seq(init_pc, end_pc, -step))))
  {
    print(pc_use)
    object1 = tryCatch(RunPCA(object = object, pc.genes = rownames(object@data), print.results = FALSE, pcs.print = NULL, do.print = F, pcs.compute = pc_use, rev.pca = T, weight.by.var = T),
                       error = function(e) {"e"}, 
                       warning = function(w) {"w"}
    )
    
    if(typeof(object1) !=  "character")
    {
      print(paste0("pc ", pc_use))
      object = object1
      return (object)
    }
  }
  return ("e")
}

compute_tsne = function(object)
{
  init_plex = max(ncol(object@raw.data)/20, 30)
  end_plex = 2
  step = (init_plex-end_plex)/10
  
  for(perplex in unique(round(seq(init_plex, end_plex, -step))))
  {
    object1 = tryCatch(RunTSNE(object = object, dims.use = 1:ncol(object@dr$pca@cell.embeddings), dim.embed = 2, check_duplicates = FALSE, perplexity = perplex),error = function(e) {"e"})
    if(typeof(object1) !=  "character")
    {
      print(paste0("perplex ", perplex))
      object = object1
      break
    }
  }
  return(object)
}

compute_clustering_min = function(object, reduction.type = "pca", gran_shresh = 2, strict_binary = T)
{
  k.param = max(round(ncol(object@raw.data)/10),3)
  k.scale = max(round(ncol(object@raw.data)/15),2)
  
  print("Clustering...")
  test_clustering_rez = function(rez)
  {
    object = FindClusters(object = object, resolution = rez, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = F)
    if(length(table(get_ident(object)))<2)
    {
      return(2)
    }
    if(length(table(get_ident(object)))>= 2)
    {
      if(min(table(get_ident(object)))<gran_shresh){
        return(NULL)
      }else if(length(table(get_ident(object))) == 2)
      {
        return(0)
      }else 
      {
        return(1)
      }
    }else{
      return(-1)
    }
  }
  
  good_rez = NULL
  direction = 1
  max_num = 3
  from = 0
  to = 30
  
  for(depth in 0:max_num)
  {
    
    if(from==to)
    {
      break
    }
    
    for(rez in seq(from, to, direction*(10^-depth)))
    {
      
      test_rez = test_clustering_rez(rez)
      print(test_rez)
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
  
  object = FindClusters(object = object, resolution = good_rez, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = F)
  return(object)
}

# compute_clustering_max = function(object, mult = 1)
# {
#   k.param = max(round(min(30, ncol(object@raw.data)/3)),3)
#   k.scale = max(round(min(25, ncol(object@raw.data)/4)),2)
#   
#   object = FindClusters(object = object, algorithm = 3, k.param = k.param, k.scale = k.scale, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T)
#   
#   if(mult<2)
#   {
#     test_clustering_rez = function(object, rez)
#     {
#         s_out = capture.output(FindClusters(object = object, resolution = rez, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T))
#       
#         s_out = unlist(strsplit(s_out[which(grepl("singletons identified", s_out))], c(" |\""), fixed = F))
#         s_out = s_out[which(grepl("^\\d+$", s_out))][1]
#     
#         if(is.null(s_out))
#         {
#           return(0)
#         }else{
#           return(as.numeric(s_out))
#         }
#     }
#     
#     rez_curr = "e"
#     
#     for(rez in seq(0.1, 10.1, 1))
#     {
#       if(test_clustering_rez(object, rez)>0)
#       {
#         rez_curr = rez
#         break
#       }
#     }
#     
#     
#     for(rez in seq(rez_curr, 0.1, -0.1))
#     {
#       if(test_clustering_rez(object, rez) == 0)
#       {
#         rez_curr = rez
#         break
#       }
#     }
#     
#     print(paste0("rez ", rez_curr))
#     object = FindClusters(object = object, resolution = rez_curr, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T)
#     object1 = ValidateClusters(object, pc.use = 1:ncol(object@dr$pca@cell.embeddings))
#     if(length(table(get_ident(object1))>1)>1)
#     {
#       print("Clustered...")
#       return (object1)
#     }else{
#         if(length(table(get_ident(object))>1)>1){
#           print("Clustered without correction...")
#           return(object)
#         }else{
#           print("singletons-allowed clustering...")
#           for(rez in seq(rez_curr, 10.1, 0.1))
#           {
#             object = FindClusters(object = object, resolution = rez, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T)
#             if(length(table(get_ident(object))>1)>1)
#             {
#               return(object)
#             }
#           }
#         }
#     }
#     return("e")
#   }else{
#     print("Unreliable clustering...")
#     for(rez in seq(0.1, 1.0, 0.1))
#     {
#       object = FindClusters(object = object, resolution = rez, k.param = k.param, k.scale = k.scale, algorithm = 3, modularity.fxn = 1, print.output = F, save.SNN = T, force.recalc = T)
#       if(length(which(table(get_ident(object))>1))>1){
#         print(paste0("rez ", rez))
#         return(object)
#       }
#     }
#   }
#   return("e")
# }

seurat_analyse_mtx = function(object, regress, regress_cc, do.magic, do.cluster = T)
{
  object = exprs(simple_preprocessing(sce_norm(object)))
  raw = NULL
  
  if(do.magic)
  {
    magic = magic_impute(as.matrix(object))
    if(!is.null(magic))
    {
      print("magic")
      raw = as(object[,colnames(magic)], "dgCMatrix")
      object = magic
    }
  }
  
  object = CreateSeuratObject(raw.data = as(2^object+1, "dgCMatrix"))
  
  return(seurat_analyse(object, regress = regress, regress_cc = regress_cc, do.cluster = do.cluster, raw_before_magic = raw))
}

seurat_analyse = function(object, regress, regress_cc, do.cluster = T, raw_before_magic = NULL)
{
  object = seur_norm_regress(object, regress = regress, regress_cc = regress_cc)
  curr_mult = get_copy_count(object)
  
  object1 = compute_pca(object)
  if(typeof(object1) !=  "character")
  {
    object = object1
    rm(object1)
  }else{
    return("e")
  }
  
  if(do.cluster)
  {
    object1 = compute_clustering_min(object)
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

seur_get_raw_and_params = function(object, ident = NULL, gran_shresh = 2, regress, regress_cc, do.magic)
{
  if(!is.null(ident) && ident  == "SeuratProject"){
    ident = NULL
  }
  #this function is required to 
  #save memory on parallel run
  #param list used in do.call
  if(is.null(ident))
  {
    raw_data = object@raw.data
  }else{
    raw_data = object@raw.data[, get_cell_names_by_ident(object, ident)]
  }
  return(list(object = raw_data, ident = ident, gran_shresh = gran_shresh, regress = regress, regress_cc = regress_cc, do.magic = do.magic))
}

seur_find_split = function(object)
{
  ident = object$ident
  gran_shresh = object$gran_shresh
  regress = object$regress
  regress_cc = object$regress_cc
  do.magic = object$do.magic
  raw_data = object$object
  
  curr_mult = get_copy_count(raw_data)
  if(
    (!is.null(ident) && (pick_last_char(toString(ident)) !=  "f")) 
    || 
    (is.null(ident) || ncol(raw_data)>= gran_shresh*2)
  )
  {
    #for(mult in seq(curr_mult, 10, 1))
    {
      #print(paste0("mult ", mult))
      object = tryCatch(seurat_analyse_mtx(raw_data, regress = regress, regress_cc = regress_cc, do.magic = do.magic, do.cluster = T),error = function(e) {print(e)
        return("e")})
      if(typeof(object) !=  "character")
      {
        if(length(which(table(get_ident(object))>gran_shresh))>1)
        {
          object = rename_clusters(object, ident)
          return(object)
        }
        if(length(table(get_ident(object))) == ncol(raw_data))
        {
          break
        }
      }
    }
  }
  
  object = CreateSeuratObject(raw.data = as(raw_data, "dgCMatrix"))
  
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


seur_stat = function(object)
{
  table(get_ident(object), get_orig_ident(object))
  
  plot_seur(object)
  
  plot_seur(object, method =  "pca")
}

split_ident = function(object, ident = NULL, gran_shresh = 2, regress = F, regress_cc = NULL, do.magic = F)
{
  param_list_for_ident = seur_get_raw_and_params(object, ident = ident, gran_shresh = gran_shresh, regress = regress, regress_cc = regress_cc, do.magic = do.magic)
  return(seur_find_split(param_list_for_ident))
}

split_all_ident = function(object, gran_shresh = 2, regress = F, regress_cc = NULL, do.magic = F)
{
  parallel_list = lapply(unique(get_ident(object)), function(ident){ seur_get_raw_and_params(object, ident = ident, gran_shresh = gran_shresh, regress = regress, regress_cc =  regress_cc, do.magic = do.magic)})
  
  #clean opened threads if the last run was terminated incorrectly
  tryCatch(doParallel::stopImplicitCluster(),error = function(e){NULL})
  
  doParallel::registerDoParallel(cores = min(parallel::detectCores()-1, length(parallel_list)))
  
  parallel_list = foreach(object = parallel_list, .options.snow = list(preschedule = TRUE)) %dopar% {
    source('Counts2Exprs.R', local = TRUE)
    object = seur_find_split(object)
    return(object)
  }
  
  doParallel::stopImplicitCluster()
  
  return(parallel_list)
}

seur_merge = function(list_seur)
{
  cl_idents = lapply(1:length(labels(list_seur)), function(i) {return(GetClusters(list_seur[[i]]))})
  new_idents = cl_idents[[1]]
  temp_object = list_seur[[1]]
  list_seur[[1]] = "NULL"
  for(i in 2:length(cl_idents))
  {
    new_idents = rbind(new_idents, cl_idents[[i]])
    temp_object = MergeSeurat(temp_object, list_seur[[i]], do.normalize = F)
    list_seur[[i]] = "NULL"
  }
  temp_object = SetClusters(temp_object, new_idents)
  temp_object = seurat_analyse(temp_object, regress = F, regress_cc = NULL, do.cluster = F)
  #temp_object = compute_tsne(temp_object)
  return (temp_object)
}





plot_seur<- function(object, method = "tsne")
{
  if(method  ==  "tsne")
  {
    object = compute_tsne(object)
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


plot_venn = function(cl_markers)
{
  VENN.LIST = cl_markers$cluster.markers.list[1:min(5, length(cl_markers$cluster.markers))]
  VENN.LIST = lapply(VENN.LIST, sort)
  
  names(VENN.LIST) = (1:min(5, length(cl_markers$cluster.markers)))-1
  grid.newpage()
  venn.plot = venn.diagram(VENN.LIST, NULL, fill = rainbow(length(names(VENN.LIST))), 
                           alpha = rep(0.5,length(names(VENN.LIST))), cex = 2, cat.fontface = 4, 
                           category.names = names(VENN.LIST), main = "DE Genes Inters")
  
  grid.draw(venn.plot)
}

find_markers = function(object, p_val = 0.05, test.use = "MAST", logfc.threshold = 0.5, only.pos = T, min.pct = 0.8, min.diff.pct = -Inf)
{
  cluster.markers = NULL
  cluster_numbers = NULL
  
  for (ident in sort(unique(get_ident(object))))
  {
    fm = tryCatch(FindMarkers(object = object, ident.1 = ident,  test.use = test.use, logfc.threshold = logfc.threshold, only.pos = only.pos, min.pct = min.pct, min.diff.pct = min.diff.pct), error = function(e){NULL})
    if(!is.null(fm))
    {
      cluster.markers = c(cluster.markers, list(ident = fm))
      cluster_numbers = c(cluster_numbers, ident)
    }
  }
  
  #filter by p_val
  for (i in 1:length(cluster.markers))
  {
    if(length(cluster.markers[i]$ident) == 0){ next }
    
    cluster.markers[i]$ident = cluster.markers[i]$ident[which(cluster.markers[i]$ident$p_val_adj<= p_val),]
  }
  
  
  cluster.markers.list = NULL
  for (i in 1:length(cluster.markers))
  {
    if(length(cluster.markers[i]$ident) == 0){
      cluster.markers.list = c(cluster.markers.list, list(item = list()))
      next }
    
    genome<-detect_genome(rownames(cluster.markers[i]$ident))
    if(!is.null(genome))
    {
      entrez_converted = gene_annot_convert_to(rownames(cluster.markers[i]$ident), org = genome$org, annot = genome$annot)
      names(entrez_converted$new) = NULL
      
      cluster.markers[i]$ident = cluster.markers[i]$ident[match(entrez_converted$old, rownames(cluster.markers[i]$ident)),]
      cluster.markers.list = c(cluster.markers.list, list(item = entrez_converted$new))
      
    }
    cluster.markers[i]$ident = cluster.markers[i]$ident[order(cluster.markers[i]$ident$avg_logFC, decreasing = T),]
  }
  cl_markers = list(cluster.markers = cluster.markers, cluster.markers.list = cluster.markers.list)
  
  return(cl_markers)
}

GO_enr = function(cl_markers, p_val = 0.05, q_val = 0.05, new_window = T, org = "m")
{
  
  cluster.markers = cl_markers$cluster.markers
  cluster.markers.list = cl_markers$cluster.markers.list
  
  names(cluster.markers.list) = rep("item", length(cluster.markers.list))
  for (i in 1:length(cluster.markers.list))
  {
    cl_mark = cluster.markers.list[i]$item
    
    x = enrichGO(gene = cl_mark, pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T, OrgDb = get_annot_db(rownames(cluster.markers$ident), org = org))
    if(nrow(x@result)>0)
    {
      x@result$Description[order(x@result$GeneRatio,decreasing = T)]
      enrichMap(x, layout = igraph::layout_components, vertex.label.cex = 1, fixed = !new_window)
      cnetplot(x, categorySize = "GeneRatio", fixed = !new_window)
    }
  }
  names(cluster.markers.list) =  sapply(1:length(cluster.markers.list), function(cl){return(paste0("X",cl-1))})
  
  res = tryCatch(compareCluster(cluster.markers.list, fun = "enrichGO", OrgDb = get_annot_db(rownames(cluster.markers$ident))), error = function(e){"e"})
  if(typeof(res) !=  "character")
  {
    plot(res, showCategory = 20)
  }
}

sort_diag = function(object)
{
  return(object[, order(sapply(1:ncol(object), function (x) {return(which.max(object[,x]))}))])
}


cluster_recursive = function(object, regress = F, regress_cc = NULL, do.magic = F)
{
  object = split_ident(CreateSeuratObject(raw.data = object, regress = regress, regress_cc = regress, do.magic = regress))
  
  while(any(sapply(sapply(get_ident(object), toString), pick_last_char)!="f"))
  {
    object = seur_merge(split_all_ident(object, regress = regress, regress_cc = regress, do.magic = regress))
  }
  return(cluster_recursive)
}




plot_seur_3d = function(object, method = "tsne", radius = 0.3, old_ident = NULL)
{
  library(plyr)
  library(scatterplot3d)
  library(rgl)
  
  object = RunTSNE(object = object, check_duplicates = FALSE, dim.embed = 3)
  
  if(method  ==  "tsne")
  {
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
  
  curr_labels = object@ident
  cells_lab = unique(curr_labels)
  labels_uniq = data.frame(names = cells_lab, colors =  rainbow(length(cells_lab)))
  
  mfrow3d(1, 2)
  next3d()
  plot3d(x = rl[,1], y = rl[,2], z = rl[,3], col = mapvalues(curr_labels, from = array(labels_uniq$names), to = array(labels_uniq$colors)), type = "s", radius = radius)
  legend3d("topright", col = array(labels_uniq$colors), legend = labels_uniq$names, cex = 2, pch = 0.1)
  next3d()
  
  if(is.null(old_ident))
  {
    curr_labels = object@meta.data$orig.ident
  }else{
    curr_labels = old_ident
  }
  
  cells_lab = unique(curr_labels)
  labels_uniq = data.frame(names = cells_lab, colors =  rainbow(length(cells_lab)))
  
  plot3d(x = rl[,1], y = rl[,2], z = rl[,3], col = mapvalues(curr_labels, from = array(labels_uniq$names), to = array(labels_uniq$colors)), type = "s", radius = radius)
  legend3d("topright", col = array(labels_uniq$colors), legend = labels_uniq$names, cex = 2, pch = 0.1)  
  
}



