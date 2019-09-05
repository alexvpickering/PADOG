#' Pathway Analysis with Down-weighting of Overlapping Genes (PADOG)
#'
#' This is a general purpose gene set analysis method that downplays the
#' importance of genes that apear often accross the sets of genes analyzed.
#'
#' The original implementation of \code{\link[PADOG]{padog}} was modified to allow two-sample
#' groups and eliminate defunct KEGG.db warning.
#'
#' @param gslist Named lists (KEGG pathway number) for each pathway, each with a character vectors of entrez ids with \code{names} attribute equal to HGNC symbols.
#' @param gs.names named character vector of KEGG pathway names with \code{names} attribute equal to KEGG pathway numbers.
#' @param rna_seq is the analysis on RNA seq data? Default is \code{FALSE}. If \code{TRUE} must supply \code{pdata}.
#' @param pdata data.frame with columns \code{lib.size} and \code{norm.factors} needed if \code{rna_seq} is \code{TRUE}.
#' @param browse Boolean set to TRUE if you want to \code{browser} through the code. For development/debugging.
#' @export
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%
#' 
#' 
#'
padog <- function (esetm = NULL, group = NULL, paired = FALSE, block = NULL, gslist = NULL, gs.names = NULL, NI = 1000, Nmin = 3, verbose = TRUE, 
                   parallel = FALSE, dseed = NULL, ncr = NULL, rna_seq = FALSE, pdata = NULL, browse = FALSE) {
  
  
  if (browse) browser()
  
  if (rna_seq) stopifnot(!is.null(pdata))
  stopifnot(is(esetm, "matrix") | is(esetm, 'Matrix'))
  # stopifnot(all(dim(esetm) > 4))
  stopifnot(is(group, "factor") | is(group, "character"))
  stopifnot(length(group) == dim(esetm)[2])
  stopifnot(all(group %in% c("c", "d")))
  # stopifnot(all(table(group) > 2))
  if (paired) {
    stopifnot(length(block) == length(group))
    stopifnot(all(table(block) == 2))
  }
  stopifnot(is(gslist, "list"))
  stopifnot(length(gslist) >= 3)
  if (!is.null(gs.names)) {
    stopifnot(length(gslist) == length(gs.names))
  }
  stopifnot(is(NI, "numeric"))
  stopifnot(NI > 5)
  
  
  stopifnot(sum(rownames(esetm) %in% as.character(unlist(gslist))) >
              10 & !any(duplicated(rownames(esetm))))
  
  Block = block
  gf = table(unlist(gslist))
  if (!all(gf == 1)) {
    if (stats::quantile(gf, 0.99) > mean(gf) + 3 * stats::sd(gf)) {
      gf[gf > stats::quantile(gf, 0.99)] <- stats::quantile(gf, 0.99)
    }
    gff <- function(x) {
      1 + ((max(x) - x)/(max(x) - min(x)))^0.5
    }
    gf = gff(gf)
    
  } else {
    fdfd = unique(unlist(gslist))
    gf = rep(1, length(fdfd))
    names(gf) <- fdfd
  }
  allGallP = unique(unlist(gslist))
  
  restg = setdiff(rownames(esetm), names(gf))
  appendd = rep(1, length(restg))
  names(appendd) <- restg
  gf = c(gf, appendd)
  stopifnot(all(!duplicated(rownames(esetm))))
  stopifnot(sum(rownames(esetm) %in% allGallP) > 10)
  if (verbose) {
    cat(paste("Starting with ", length(gslist), " gene sets!", sep = ""))
    cat("\n")
  }
  
  # drop pathways with less than Nmin genes
  gslist = gslist[unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x)) >= Nmin
  }))]
  gs.names = gs.names[names(gslist)]
  stopifnot(length(gslist) >= 3)
  if (verbose) {
    cat(paste("Analyzing ", length(gslist), " gene sets with ", Nmin, " or more genes!", sep = ""))
    cat("\n")
  }
  if (!is.null(dseed))
    set.seed(dseed)
  G = factor(group)
  Glen = length(G)
  tab = table(G)
  idx = which.min(tab)
  minG = names(tab)[idx]
  minGSZ = tab[idx]
  bigG = rep(setdiff(levels(G), minG), length(G))
  block = factor(Block)
  topSigNum = dim(esetm)[1]
  combFun = function(gi, countn = TRUE) {
    g = G[gi]
    tab = table(g)
    if (countn) {
      minsz = min(tab)
      ifelse(minsz > 10, -1, choose(length(g), minsz))
      
    } else {
      dup = which(g == minG)
      cms = utils::combn(length(g), tab[minG])
      del = apply(cms, 2, setequal, dup)
      if (paired) {
        cms = cms[, order(del, decreasing = TRUE), drop = FALSE]
        cms[] = gi[c(cms)]
        cms
        
      } else {
        cms[, !del, drop = FALSE]
      }
    }
  }
  if (paired) {
    bct = tapply(seq_along(G), block, combFun, simplify = TRUE)
    nperm = ifelse(any(bct < 0), -1, prod(bct))
    if (nperm < 0 || nperm > NI) {
      btab = tapply(seq_along(G), block, `[`, simplify = FALSE)
      bSamp = function(gi) {
        g = G[gi]
        tab = table(g)
        bsz = length(g)
        minsz = tab[minG]
        cms = do.call(cbind, replicate(NI, sample.int(bsz, minsz), simplify = FALSE))
        cms[] = gi[c(cms)]
        cms
      }
      combidx = do.call(rbind, lapply(btab, bSamp))
      
    } else {
      bcomb = tapply(seq_along(G), block, combFun, countn = FALSE, simplify = FALSE)
      colb = expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1, , drop = FALSE]
      combidx = mapply(function(x, y) x[, y, drop = FALSE], bcomb, colb, SIMPLIFY = FALSE)
      combidx = do.call(rbind, combidx)
    }
  } else {
    nperm = combFun(seq_along(G))
    if (nperm < 0 || nperm > NI) {
      combidx = do.call(cbind, replicate(NI, sample.int(Glen, minGSZ), simplify = FALSE))
    } else {
      combidx = combFun(seq_along(G), countn = FALSE)
    }
  }
  NI = ncol(combidx)
  
  deINgs = intersect(rownames(esetm), unlist(gslist))
  gslistINesetm = lapply(gslist, match, table = deINgs, nomatch = 0)
  MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI + 1)
  gsScoreFun <- function(G, block) {
    # these two arguments are needed in parallel computing for the environment in
    # model.matrix
    force(G)
    force(block)
    if (ite > 1) {
      G = bigG
      G[combidx[, ite - 1]] = minG
      G = factor(G)
    }
    if (paired) {
      design <- stats::model.matrix(~0 + G + block)
      colnames(design) <- substr(colnames(design), 2, 100)
      
    } else {
      design <- stats::model.matrix(~0 + G)
      colnames(design) <- levels(G)
    }
    
    fit2 <- fit_ebayes(esetm, contrasts = 'd-c', mod = design, pdata = pdata, rna_seq = rna_seq)
    
    aT1 <- limma::topTable(fit2, coef = 1, number = topSigNum)
    aT1$ID = rownames(aT1)
    de = abs(aT1$t)
    names(de) <- aT1$ID
    degf = scale(cbind(de, de * gf[names(de)]))
    rownames(degf) = names(de)
    degf = degf[deINgs, , drop = FALSE]
    sapply(gslistINesetm, function(z) {
      X = stats::na.omit(degf[z, , drop = FALSE])
      colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
    })
  }
  if (parallel && requireNamespace("doParallel", quietly = TRUE) &&
      requireNamespace("parallel", quietly = TRUE)) {
    ncores = parallel::detectCores()
    if (!is.null(ncr))
      ncores = min(ncores, ncr)
    
    clust = parallel::makeCluster(ncores)
    
    doParallel::registerDoParallel(clust)
    tryCatch({
      parRes = foreach::foreach(ite = 1:(NI + 1), .combine = "c",
                                .packages = "limma") %dorng% {
                                  Sres <- gsScoreFun(G, block)
                                  tmp <- list(t(Sres))
                                  names(tmp) <- ite
                                  if (verbose && (ite%%10 == 0)) {
                                    cat(ite, "/", NI, "\n")
                                  }
                                  tmp
                                }
      parRes = do.call(cbind, parRes[order(as.numeric(names(parRes)))])
      evenCol = (1:ncol(parRes))%%2 == 0
      MSabsT[] = parRes[, !evenCol]
      MSTop[] = parRes[, evenCol]
      rm(parRes)
    }, finally = parallel::stopCluster(clust))
    
  } else {
    if (parallel)
      message("Execute in serial! Packages 'doParallel' and 'parallel' \n needed for parallelization!")
    for (ite in 1:(NI + 1)) {
      Sres <- gsScoreFun(G, block)
      MSabsT[, ite] <- Sres[1, ]
      MSTop[, ite] <- Sres[2, ]
      if (verbose && (ite %% 10 == 0)) {
        cat(ite, "/", NI, "\n")
      }
    }
  }
  meanAbsT0 = MSabsT[, 1]
  padog0 = MSTop[, 1]
  plotIte = min(NI, 21)
  MSabsT_raw = MSabsT
  MSTop_raw = MSTop
  
  # standardize scores
  MSabsT = scale(MSabsT)
  MSTop = scale(MSTop)
  
  # compute p-values
  mff = function(x) {
    mean(x[-1] > x[1], na.rm = TRUE)
  }
  PSabsT = apply(MSabsT, 1, mff)
  PSTop = apply(MSTop, 1, mff)
  PSabsT[PSabsT == 0] <- 1/NI/100
  PSTop[PSTop == 0] <- 1/NI/100
  
  # estimate FDR
  pval = list()
  fdrs = lapply(list(MSabsT, MSTop), function(x) {
    p1 = rowMeans(x[,-1,drop=FALSE] > x[,1], na.rm=TRUE)
    x = x[,-1,drop=FALSE]
    p0 = sapply(1:ncol(x), function(z) rowMeans(x[,-z,drop=FALSE] > x[,z], na.rm=TRUE))
    pval <<- c(pval, list(p0))
    getFDR(p0, p1)
  })
  names(fdrs) = c("FDRmeanAbsT", "FDRpadog")
  names(pval) = c("AbsmT", "PADOG") #use method names
  
  if (!is.null(gs.names)) {
    myn = gs.names
    
  } else {
    myn = names(gslist)
  }
  SIZE = unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x))
  }))
  res = data.frame(Name = myn, ID = names(gslist), Size = SIZE,
                   meanAbsT0, padog0, PmeanAbsT = PSabsT, fdrs, Ppadog = PSTop,
                   stringsAsFactors = FALSE)
  ord = order(res$Ppadog, -res$padog0)
  res = res[ord, ]
  pval = lapply(pval, function(x) {x = x[ord,,drop=FALSE]; rownames(x) = res$ID; x})
  list(res=res, pval=pval)
}


#' Perform eBayes analysis from limma.
#'
#' Generates contrast matrix then runs eBayes analysis from limma.
#' 
#' @inheritParams padog 
#' @param contrasts comparison to make. For \code{\link{padog}} is \code{'d-c'}.
#' @param mod design matrix
#'
#' @return result from call to limma \code{eBayes}.
#' @export
#' @keywords internal
fit_ebayes <- function(esetm, contrasts, mod, pdata, rna_seq = FALSE) {
  
  if (rna_seq) {
    lib.size <- pdata$lib.size * pdata$norm.factors
    v <- limma::voom(esetm, mod, lib.size)
    fit  <- limma::lmFit(v)
    
  } else {
    fit <- limma::lmFit(esetm, mod)
  }
  
  cont.matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  return(limma::eBayes(fit2))
}