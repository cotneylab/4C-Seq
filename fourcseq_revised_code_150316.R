library(GenomicAlignments)

countFragmentOverlaps_edi<-function (object, trim = 0, minMapq = 0, shift = 0, overlaptype="any") 
{
    stopifnot(class(object) == "FourC")
    if (length(rowRanges(object)) == 0) 
        stop("Add fragments before calling 'findViewpointFragments'")
    cat("reading bam files\n")
    bamFiles = file.path(metadata(object)$bamFilePath, colData(object)$bamFile)
    colData(object)$orignialReads = sapply(bamFiles, function(bamFile) countBam(bamFile)$records)
    reads = lapply(bamFiles, function(bamfile) {
        what <- c("mapq")
        flag <- scanBamFlag(isUnmappedQuery = FALSE)
        param <- ScanBamParam(what = what, flag = flag)
        readGAlignments(bamfile, param = param)
    })
    colData(object)$rawReads = sapply(reads, length)
    if (minMapq >= 0) 
        reads = lapply(reads, function(rga, minMapq) rga[mcols(rga)$mapq > 
            minMapq], minMapq = minMapq)
    colData(object)$lowQualityReads = colData(object)$rawReads - 
        sapply(reads, length)
    reads = lapply(reads, granges)
    if (trim > 0) {
        reads = lapply(reads, function(gr, trim) {
            start(gr[strand(gr) == "+"]) = start(gr[strand(gr) == 
                "+"]) + trim
            end(gr[strand(gr) == "-"]) = end(gr[strand(gr) == 
                "-"]) - trim
            gr
        }, trim = trim)
    }
    cat("calculating overlaps\n")
    frag <- rowRanges(object)
    
    if(overlaptype == "any") {
		strand(frag) <- "+"
		countsLeftFragmentEnd <- sapply(reads, countOverlaps, query = frag, 
			type = overlaptype, maxgap = shift,ignore.strand =F) #JY edited here, to count reads from both cutting sites
		strand(frag) <- "-"
		countsRightFragmentEnd <- sapply(reads, countOverlaps, query = frag, 
			type = overlaptype , maxgap = shift,ignore.strand =F) #JY edited here, to count reads from both cutting sites
	} else {
		#go back to the default method #JY
		strand(frag) <- "+"
		countsLeftFragmentEnd <- sapply(reads, countOverlaps, query = frag, 
			type = c("start"), maxgap = shift)
		strand(frag) <- "-"
		countsRightFragmentEnd <- sapply(reads, countOverlaps, query = frag, 
			type = c("end"), maxgap = shift)
	}      
    
    colData(object)$mappedReads = apply(countsLeftFragmentEnd, 
        2, sum) + apply(countsRightFragmentEnd, 2, sum)
    colData(object)$mappingRatio = colData(object)$mappedReads/(colData(object)$rawReads - 
        colData(object)$lowQualityReads)
    metaDataFrame <- DataFrame(type = rep("countInfo", 5), description = rep("", 
        5))
    idx <- colnames(colData(object)) %in% c("orignialReads", 
        "rawReads", "lowQualityReads", "mappedReads", "mappingRatio")
    mcols(colData(object))[idx, ] <- metaDataFrame
    assays(object) <- SimpleList(countsLeftFragmentEnd = countsLeftFragmentEnd, 
        countsRightFragmentEnd = countsRightFragmentEnd)
    assays(object) <- SimpleList(countsLeftFragmentEnd = countsLeftFragmentEnd, 
        countsRightFragmentEnd = countsRightFragmentEnd)
    object
}

getZScores_edi <-function (object, removeZeros = TRUE, minCount = 40, maxDist=NULL, minDist = NULL, 
    fitFun = "distFitMonotoneSymmetric", sdFun = mad, ...) 
{
    stopifnot(class(object) == "FourC")
    if (!c("counts") %in% names(assays(object))) 
        stop("No assay 'counts' found. Use combineFragEnds first.")
    if (c("zScore") %in% names(assays(object))) 
        stop("z-scores are already calculated. To recalculate z-scores use the object returned by 'combineFragEnds'.")
    viewpoint = unique(colData(object)$viewpoint)
    if (length(viewpoint) != 1) 
        stop("None or more than one viewpoint are contained in the 'FourC' object.\n         Use a 'FourC' object that contains only one viewpoint.")
    print(viewpoint)
    object <- object[seqnames(object) == unique(colData(object)$chr), 
        ]
    fragData = getDistAroundVp(viewpoint, colData(object), rowRanges(object))
    medianCounts <- apply(counts(object), 1, median)
    if (!is.null(minDist)) {
        tooClose = abs(fragData$dist) <= minDist
    }
    else {
        toLeft <- fragData$dist > -20000 & fragData$dist < 0 & 
            !is.na(fragData$dist)
        afterMin <- 1:sum(toLeft) > tail(which(sign(diff(medianCounts[toLeft])) < 
            0), 1) + 1
        toExclude <- which(toLeft)[afterMin]
        toRight <- fragData$dist < 20000 & fragData$dist > 0 & 
            !is.na(fragData$dist)
        beforeMin <- 1:sum(toRight) < which(sign(diff(medianCounts[toRight])) > 
            0)[1]
        toExclude <- c(toExclude, which(toRight)[beforeMin])
        toExclude = union(toExclude, which(abs(fragData$dist) < 
            1000))
        tooClose = rep(FALSE, length(fragData$dist))
        tooClose[toExclude] = TRUE
    }
    
    #edited by JY, 03122015
    if(!is.null(maxDist)) {
    	tooFar = (abs(fragData$dist) > maxDist)
    } else {
    	tooFar= rep(0,length(abs(fragData$dist)))
    }
    
    #remove fragment too close to the viewpoint
    fragData$posLeft[tooClose] = FALSE
    fragData$posRight[tooClose] = FALSE
    fragData$tooClose = tooClose
    
    #remove fragments too far away from the viewpoint #edited by JY, 03122015
    fragData$posLeft[tooFar] = FALSE
    fragData$posRight[tooFar] = FALSE
	fragData$tooFar = tooFar
    
    lowCounts = medianCounts < minCount
    fragData$posLeft[lowCounts] = FALSE
    fragData$posRight[lowCounts] = FALSE
    fragData$lowCounts = lowCounts
    fragData$selectedForFit <- (fragData$posLeft | fragData$posRight)
    newCols <- c("tooClose", "tooFar", "lowCounts", "selectedForFit")
    mcolsRows <- DataFrame(type = rep("fragmentSelection", length(newCols)), 
        description = rep("", length(newCols)))
    mcols(mcols(fragData))[colnames(mcols(fragData)) %in% newCols, 
        ] <- mcolsRows
    rowRanges(object) <- fragData
    fragData = as.data.frame(fragData)
    colData(object)$condition <- factor(colData(object)$condition)
    dds <- object[mcols(object)$selectedForFit, ]
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    if (attr(dispersionFunction(dds), "fitType") != "parametric") {
        stop("Failed to estimate the parameters of the Variance stabilizing transformation.")
    }
    else {
        coefs <- attr(dispersionFunction(dds), "coefficients")
        attr(dds, "vst") <- function(q) {
            log((1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * 
                q + 2 * sqrt(coefs["asymptDisp"] * q * (1 + coefs["extraPois"] + 
                coefs["asymptDisp"] * q)))/(4 * coefs["asymptDisp"]))/log(2)
        }
        attr(dds, "inverse-vst") <- function(q) {
            (4 * coefs["asymptDisp"] * 2^q - (1 + coefs["extraPois"]))^2/(4 * 
                coefs["asymptDisp"] * (1 + coefs["extraPois"] + 
                (4 * coefs["asymptDisp"] * 2^q - (1 + coefs["extraPois"]))))
        }
    }
    trafo <- FourCSeq:::getVST(dds)(counts(dds))
    fit <- apply(trafo, 2, fitFun, fragData = as.data.frame(rowRanges(dds)), 
        removeZeros = removeZeros, ...)
    residuals <- trafo - fit
    sd <- apply(residuals, 2, sdFun)
    zScore <- sweep(residuals, 2, sd, "/")
    pValue <- apply(zScore, 2, pnorm, lower.tail = FALSE)
    pAdjusted <- apply(pValue, 2, p.adjust, method = "BH")
    colData(dds)$sd <- sd
    idx <- which(colnames(colData(dds)) == "sd")
    metaDataFrame <- DataFrame(type = "intermediate", description = "sd/mad calculated from the residuals")
    mcols(colData(dds))[idx, ] <- metaDataFrame
    metadata(dds)$parameter <- DataFrame(fitFun = fitFun, removeZeros = removeZeros, 
        minCount = minCount, sdFun = sdFun, ...)
    if (!is.null(minDist)) 
        metadata(dds)$parameter$minDist = minDist
    assays(dds) <- c(assays(dds), SimpleList(trafo = trafo, fit = fit, 
        zScore = zScore, pValue = pValue, pAdjusted = pAdjusted))
    invisible(dds)
}

