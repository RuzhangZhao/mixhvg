
comb_rank<-function(input_lst){
  input_lst_order<-list()
  for(i in 1:length(input_lst)){
    input_lst_order[[i]]<-order(order(input_lst[[i]],decreasing = FALSE))
  }

  apply(matrix(unlist(input_lst_order),
               ncol=length(input_lst_order),byrow = FALSE),1,FUN = max)
}

FindFeatureVal<-function(method.names,
                         counts=NULL,
                         normalizedcounts=NULL,
                         lognormalizedcounts=NULL,
                         PFlog1pPF=NULL,
                         loess.span = 0.3,
                         clip.max = "auto",
                         num.bin = 20,
                         binning.method = "equal_width",
                         verbose = FALSE){
  switch(method.names,
         "seuratv3"={
           feature_val <- FindVariableFeatures(counts,
                                               loess.span = loess.span,
                                               clip.max = clip.max,
                                               num.bin = num.bin,
                                               binning.method = binning.method,
                                               verbose = verbose)$vst.variance.standardized
         },
         "seuratv1"={
           feature_val <- FindVariableFeatures(
             lognormalizedcounts,
             selection.method = "disp",
             num.bin = num.bin,
             binning.method = binning.method,
             verbose = verbose)$mvp.dispersion
           feature_val[is.na(feature_val)]<-0
           feature_val[feature_val<0]<-0
         },
         "mv_ct"={
           sce <- SingleCellExperiment(list(counts=counts))
           sce@assays@data$logcounts<-sce@assays@data$counts
           dec <- modelGeneVar(sce)
           dec.var <- dec@listData$bio
           dec.keep <- !is.na(dec.var) & dec.var > 0
           dec.var[!dec.keep]<-0
           feature_val<-dec.var
         },
         "mv_nc"={
           sce <- SingleCellExperiment(list(counts=normalizedcounts))
           sce@assays@data$logcounts<-sce@assays@data$counts
           dec <- modelGeneVar(sce)
           dec.var <- dec@listData$bio
           dec.keep <- !is.na(dec.var) & dec.var > 0
           dec.var[!dec.keep]<-0
           feature_val<-dec.var
         },
         "mv_PFlogPF"={
           sce <- SingleCellExperiment(list(counts=PFlog1pPF))
           sce@assays@data$logcounts<-sce@assays@data$counts
           dec <- modelGeneVar(sce)
           dec.var <- dec@listData$bio
           dec.keep <- !is.na(dec.var) & dec.var > 0
           dec.var[!dec.keep]<-0
           feature_val<-dec.var
         },
         "scran"={
           sce <- SingleCellExperiment(list(counts=counts))
           sce <- logNormCounts(sce)
           dec <- modelGeneVar(sce)
           dec.var <- dec@listData$bio
           dec.keep <- !is.na(dec.var) & dec.var > 0
           dec.var[!dec.keep]<-0
           feature_val<-dec.var
         },
         "scran_pos"={
           sce <- SingleCellExperiment(list(counts=counts))
           sce <- logNormCounts(sce)
           dec<- modelGeneVarByPoisson(sce)
           dec.var <- dec@listData$bio
           dec.keep <- !is.na(dec.var) & dec.var > 0
           dec.var[!dec.keep]<-0
           feature_val<-dec.var
         },
         "logmv_nc"={
           feature_val <- FindVariableFeatures(normalizedcounts,
                                               loess.span = loess.span,
                                               clip.max = clip.max,
                                               num.bin = num.bin,
                                               binning.method = binning.method,
                                               verbose = verbose)$vst.variance.standardized
         },
         "logmv_lognc"={
           feature_val <- FindVariableFeatures(lognormalizedcounts,
                                               loess.span = loess.span,
                                               clip.max = clip.max,
                                               num.bin = num.bin,
                                               binning.method = binning.method,
                                               verbose = verbose)$vst.variance.standardized
         },
         "logmv_PFlogPF"={
           feature_val <- FindVariableFeatures(PFlog1pPF,
                                               loess.span = loess.span,
                                               clip.max = clip.max,
                                               num.bin = num.bin,
                                               binning.method = binning.method,
                                               verbose = verbose)$vst.variance.standardized
         },
         "disp_PFlogPF"={
           feature_val <- FindVariableFeatures(
             PFlog1pPF,
             selection.method = "disp",
             num.bin = num.bin,
             binning.method = binning.method,
             verbose = verbose)$mvp.dispersion
           feature_val[is.na(feature_val)]<-0
           feature_val[feature_val<0]<-0
         },
         "mean_max_ct"={
           feature_val<-rowMeans(counts)
         },
         "mean_max_nc"={
           feature_val<-rowMeans(normalizedcounts)
         },
         "mean_max_lognc"={
           feature_val<-rowMeans(lognormalizedcounts)
         },
         {
           print("wrong input!")
         }
  )
  return(feature_val)
}

#' FindVariableFeaturesMix
#'
#' @details The function inherits from FindVariableFeatures function of Seurat Package. Refer to \url{https://github.com/RuzhangZhao/mixhvg} for user manual.
#'
#' @param object An object, SeuratObject and matrix(including sparse matrix) are both acceptable
#' @param method.names The following methods can be directly used for highly variable feature selection. The mixture of methods take a vector of method list, e.g. c("mv_nc","scran_pos","seuratv1"), which is also default.
#' \itemize{
#' \item{scran: }{Use mean-variance curve adjustment on lognormalized count matrix, which is scran ModelGeneVar.}
#' \item{mv_ct: }{Use mean-variance curve adjustment on count matrix, inherited from scran ModelGeneVar.}
#' \item{mv_nc: }{Use mean-variance curve adjustment on normalized count matrix, inherited from scran ModelGeneVar.}
#' \item{mv_lognc: }{The same as scran.}
#' \item{mv_PFlogPF: }{Use mean-variance curve adjustment on PFlog1pPF matrix, inherited from scran ModelGeneVar.}
#' \item{scran_pos: }{Use scran poisson version, modelGeneVarByPoisson.}
#' \item{seuratv3: }{Use logmean-logvariance curve adjustment on count matrix, which is vst, Seurat FindVariableFeatures Function(\url{https://satijalab.org/seurat/reference/findvariablefeatures}).}
#' \item{logmv_ct: }{The same as seuratv3.}
#' \item{logmv_nc: }{Use logmean-logvariance curve adjustment on normalized count matrix, inherited from seuratv3(vst).}
#' \item{logmv_lognc: }{Use logmean-logvariance curve adjustment on lognormalized count matrix, inherited from seuratv3(vst).}
#' \item{logmv_PFlogPF: }{Use logmean-logvariance curve adjustment on PFlog1pPF matrix, inherited from seuratv3(vst).}
#' \item{seuratv1: }{Use dispersion on lognormalized count matrix, which is dispersion (disp), Seurat FindVariableFeatures Function(\url{https://satijalab.org/seurat/reference/findvariablefeatures}).}
#' \item{disp_lognc: }{The same as seuratv1.}
#' \item{disp_PFlogPF: }{Use dispersion on PFlog1pPF matrix, inherited from seuratv1(disp).}
#' \item{mean_max_ct: }{Highly Expressed Features with respect to count matrix.}
#' \item{mean_max_nc: }{Highly Expressed Features with respect to normalized count matrix.}
#' \item{mean_max_lognc: }{Highly Expressed Features with respect to lognormalized count matrix}
#' }
#' @param nfeatures Number of features to select as top variable features.
#' @param loess.span (Only work for logmv based methods like seuratv3). Loess span parameter used when fitting the variance-mean relationship
#' @param clip.max (Only work for logmv based methods like seuratv3). After standardization values larger than clip.max will be set to clip.max; default is 'auto' which sets this value to the square root of the number of cells
#' @param num.bin (Only work for logmv or dispersion based methods)Total number of bins to use in the scaled analysis (default is 20)
#' @param binning.method Specifies how the bins should be computed. Available methods are:
#' \itemize{
#' \item{equal_width: }{each bin is of equal width along the x-axis[default].}
#' \item{equal_frequency: }{each bin contains an equal number of features (can increase statistical power to detect overdispersed features at high expression values, at the cost of reduced resolution along the x-axis).}
#' }
#' @param verbose Whether to show progress bar for calculations. Default is FALSE.
#'
#'
#'
#' @return object:  If the input is SeuratObject, the return is also SeuratObject; if the input is matrix(including sparse matrix), the return is the highly variable feature names.
#'
#'
#' @import Matrix
#' @import scran
#' @import Seurat
#' @importFrom Seurat NormalizeData FindVariableFeatures VariableFeatures DefaultAssay
#' @importFrom scran modelGeneVar modelGeneVarByPoisson
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @importFrom methods as
#' @export
#'
#' @examples
#' if(0){
#' simple_matrix<-matrix(1:2e4,nrow=4000,ncol=5)
#' rownames(simple_matrix)<-1:nrow(simple_matrix)
#' colnames(simple_matrix)<-1:ncol(simple_matrix)
#' simple_matrix_HVG<-FindVariableFeaturesMix(simple_matrix)
#' }
#'
FindVariableFeaturesMix<-function(object,
                                method.names = c("mv_nc","scran_pos","seuratv1"),
                                nfeatures = 2000,
                                loess.span = 0.3,
                                clip.max = "auto",
                                num.bin = 20,
                                binning.method = "equal_width",
                                verbose = FALSE){
  if (nrow(object) < nfeatures){
    stop("nfeatures should be smaller than
      the number of features in expression
      matrix")
  }
  if(is.null(rownames(object)[1])){
    rownames(object)<-c(1:nrow(object))
  }
  if(is.null(colnames(object)[1])){
    colnames(object)<-c(1:ncol(object))
  }else if(length(unique(colnames(object)))<
           ncol(object) ) {
    print("WARN: There are duplicated cell names! Make cell names unique by renaming!")
    colnames(object)<-make.unique(colnames(object))
  }else if(length(unique(rownames(object)))<
           nrow(object) ) {
    print("WARN: There are duplicated gene names! Make gene names unique by renaming!")
    rownames(object)<-make.unique(rownames(object))
  }

  if(inherits(x = object, 'Seurat')){
    res_return<-"Return Object"
    counts<-object@assays[[DefaultAssay(object)]]@counts
  }else if(inherits(x = object, 'Matrix') | inherits(x = object, 'matrix')){
    if (!inherits(x = object, what = 'dgCMatrix')) {
      object <- as(object = object, Class = 'dgCMatrix')
    }
    res_return<-"Return Features"
    counts<-object
  }else{
    stop("Input only accept SeuratObject or matrix(including sparse)!")
  }
  method.names[method.names == "disp_lognc"]<-"seuratv1"
  method.names[method.names == "logmv_ct"]<-"seuratv3"
  method.names[method.names == "mv_lognc"]<-"scran"
  method.names<-unique(method.names)
  pf_group<-c("disp_PFlogPF","logmv_PFlogPF","mv_PFlogPF")
  ln_group<-c("mean_max_lognc","logmv_lognc","seuratv1")
  nc_group<-c("mean_max_nc","logmv_nc","mv_nc")
  ct_group<-c("mean_max_ct","seuratv3","mv_ct","scran_pos","scran")
  normalizedcounts<-NULL
  lognormalizedcounts<-NULL
  PFlog1pPF<-NULL
  if(sum(method.names%in%nc_group)>0){
    if (!inherits(x = counts, 'Matrix')) {
      counts <- as(object = as.matrix(x = counts), Class = 'Matrix')
    }
    if (!inherits(x = counts, what = 'dgCMatrix')) {
      counts <- as(object = counts, Class = 'dgCMatrix')
    }
    lognormalizedcounts<-NormalizeData(counts,verbose=FALSE)
    normalizedcounts<-lognormalizedcounts
    normalizedcounts@x<-exp(normalizedcounts@x)-1
    #lognormalizedcounts<-as.matrix(lognormalizedcounts)
    #normalizedcounts<-as.matrix(normalizedcounts)
    #counts<-as.matrix(counts)
  }else if(sum(method.names%in%ln_group)>0){
    lognormalizedcounts<-NormalizeData(counts,verbose=FALSE)
    #lognormalizedcounts<-as.matrix(lognormalizedcounts)
  }
  if(sum(method.names%in%pf_group)>0){
    PFlog1pPF<-t(t(counts)/colSums(counts))*mean(colSums(counts))
    PFlog1pPF<-log1p(PFlog1pPF)
    PFlog1pPF<-t(t(PFlog1pPF)/colSums(PFlog1pPF))*mean(colSums(PFlog1pPF))
    #PFlog1pPF<-as.matrix(PFlog1pPF)
  }

  if(length(method.names) == 1){
    feature_val<-FindFeatureVal(method.names,
                                counts=counts,
                                normalizedcounts = normalizedcounts,
                                lognormalizedcounts = lognormalizedcounts,
                                PFlog1pPF = PFlog1pPF,
                                loess.span = loess.span,
                                clip.max = clip.max,
                                num.bin = num.bin,
                                binning.method = binning.method,
                                verbose = verbose)
  }else{
    feature_val_list<-list()
    for(i in 1:length(method.names)){
      feature_val_list[[i]]<-FindFeatureVal(method.names[i],
                                            counts=counts,
                                            normalizedcounts = normalizedcounts,
                                            lognormalizedcounts = lognormalizedcounts,
                                            PFlog1pPF = PFlog1pPF,
                                            loess.span = loess.span,
                                            clip.max = clip.max,
                                            num.bin = num.bin,
                                            binning.method = binning.method,
                                            verbose = verbose)
      feature_val<-comb_rank(feature_val_list)
    }
  }
  if (res_return == "Return Object"){
    VariableFeatures(object)<-rownames(counts)[order(feature_val,decreasing = TRUE)[1:nfeatures]]
    return(object)
  }
  if (res_return == "Return Features"){
    return(rownames(counts)[order(feature_val,decreasing = TRUE)[1:nfeatures]])
  }
}
