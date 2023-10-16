SelectIntegrationFeaturesMix <- function(
        object.list,
        nfeatures = 2000,
        assay = NULL,
        verbose = TRUE,
        fvf.nfeatures = 2000,
        ...
) {
    if (!is.null(x = assay)) {
        if (length(x = assay) != length(x = object.list)) {
            stop("If specifying the assay, please specify one assay per object in the object.list")
        }
        for (ii in length(x = object.list)) {
            DefaultAssay(object = object.list[[ii]]) <- assay[ii]
        }
    } else {
        assay <- sapply(X = object.list, FUN = DefaultAssay)
    }
    for (ii in 1:length(x = object.list)) {
        if (length(x = VariableFeatures(object = object.list[[ii]])) == 0) {
            if (verbose) {
                message(paste0("No variable features found for object", ii, " in the object.list. Running FindVariableFeatures ..."))
            }
            object.list[[ii]] <- FindVariableFeaturesMix(object = object.list[[ii]], nfeatures = fvf.nfeatures, verbose = verbose, ...)
        }
    }
    var.features <- unname(obj = unlist(x = lapply(
        X = 1:length(x = object.list),
        FUN = function(x) VariableFeatures(object = object.list[[x]], assay = assay[x]))
    ))
    var.features <- sort(x = table(var.features), decreasing = TRUE)
    for (i in 1:length(x = object.list)) {
        var.features <- var.features[names(x = var.features) %in% rownames(x = object.list[[i]][[assay[i]]])]
    }
    tie.val <- var.features[min(nfeatures, length(x = var.features))]
    features <- names(x = var.features[which(x = var.features > tie.val)])
    vf.list <- lapply(X = object.list, FUN = VariableFeatures)
    if (length(x = features) > 0) {
        feature.ranks <- sapply(X = features, FUN = function(x) {
            ranks <- sapply(X = vf.list, FUN = function(vf) {
                if (x %in% vf) {
                    return(which(x = x == vf))
                }
                return(NULL)
            })
            median(x = unlist(x = ranks))
        })
        features <- names(x = sort(x = feature.ranks))
    }
    features.tie <- var.features[which(x = var.features == tie.val)]
    tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
        ranks <- sapply(X = vf.list, FUN = function(vf) {
            if (x %in% vf) {
                return(which(x = x == vf))
            }
            return(NULL)
        })
        median(x = unlist(x = ranks))
    })
    features <- c(
        features,
        names(x = head(x = sort(x = tie.ranks), nfeatures - length(x = features)))
    )
    return(features)
}
