#' @useDynLib velocyto.R
#' @import MASS
#' @import stats
#' @import graphics
#' @importFrom Matrix colSums rowSums spMatrix Diagonal t rowMeans colMeans rowSums colSums diag
#' @importFrom utils read.delim
#' @importFrom pcaMethods pca
#' @importFrom mgcv gam s
#' @importFrom parallel mclapply
#' @importFrom cluster pam
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices adjustcolor colorRampPalette
#//' @importFrom Biostrings PWM matchPWM
#//' @importFrom GenomicRanges GRanges
#//' @importFrom IRanges IRanges
#//' @importFrom data.table fread
#' @importFrom methods as
NULL

# optional imports
# @import igraph
# @importFrom abind abind
# @import hdf5r
# @importFrom edgeR calcNormFactors
# @import GenomicAlignments
# @import Rsamtools
# @importFrom Rtsne Rtsne


##' Estimate RNA velocity using gene-relative slopes
##'
##' @param emat - spliced (exonic) count matrix
##' @param nmat - unspliced (nascent) count matrix
##' @param deltaT - amount of time to project the cell forward
##' @param smat - optional spanning read matrix (used in offset calculations)
##' @param steady.state.cells - optional set of steady-state cells on which the gamma should be estimated (defaults to all cells)
##' @param kCells - number of k nearest neighbors (NN) to use in slope calculation smoothing
##' @param cellKNN - optional pre-calculated cell KNN matrix
##' @param kGenes - number of genes (k) to use in gene kNN pooling
##' @param old.fit - optional old result (in this case the slopes and offsets won't be recalculated, and the same kNN graphs will be used)
##' @param mult - library scaling factor (1e6 in case of FPM)
##' @param min.nmat.smat.correlation - minimum required Spearman rank correlation between n and s counts of a gene
##' @param min.nmat.emat.correlation - minimum required Spearman rank correlation between n and e counts of a gene
##' @param min.nmat.emat.slope - minimum sloope of n~e regression
##' @param zero.offset - should offset be set to zero, or determined (through smat regression or using near-0 e cases)
##' @param deltaT2 - scaling of the projected difference vector (normally should be set to 1)
##' @param fit.quantile perform gamma fit on a top/bottom quantiles of expression magnitudes
##' @param diagonal.quantiles whether extreme quantiles should be computed diagonally
##' @param show.gene an optional name of a gene for which the velocity estimation details should be shown (instead of estimating all velocities)
##' @param do.par whether the graphical device parameters should be reset as part of show.gene (default=TRUE)
##' @param cell.dist - cell distance to use in cell kNN pooling calculations
##' @param emat.size - pre-calculated cell sizes for the emat (spliced) matrix
##' @param nmat.size - pre-calculated cell sizes for the nmat (unspliced) matrix
##' @param cell.emb - cell embedding to be used in show.gene function
##' @param cell.colors - cell colors to be used in show.gene function
##' @param expression.gradient - color palette used to show the expression magnitudes in show.gene function
##' @param residual.gradient - color palette used to show the u residuals in show.gene function
##' @param n.cores - number of cores to use
##' @param verbose - output messages about progress
##' @return a list with velocity results, including the current normalized expression state ($current), projected ($projected) over a certain time ($deltaT), unscaled transcriptional change ($deltaE), fit results ($gamma, $ko, $sfit if spanning reads were used), optional cell pooling parameters ($cellKNN, $kCells), kNN-convolved normalized matrices (conv.nmat.norm and conv.emat.norm), library scale ($mult)
##' @examples
##' \dontrun{
##'  # use min/max quantile gamma fit (recommended option when one can afford to do cell kNN smoothing)
##'  # The example below uses k=5 cell kNN pooling, and top/bottom 2% exprssion quantiles
##'  # emat and nmat are spliced (exonic) and unspliced (intronic) molecule/read count matirces
##' (preferably filtered for informative genes)
##'  rvel <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02)
##'
##'  # alternativly, the function can be used to visualize gamma fit and regression for a
##' particular gene. here we pass embedding (a matrix/data frame with rows named with cell names,
##' and columns corresponding to the x/y coordinates)
##' 
##'  # and cell colors. old.fit is used to save calculation time.
##'  gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02,
##'     old.fit=rvel,show.gene='Chga',cell.emb=emb,cell.colors=cell.colors)
##' }
##' @export
gene.relative.velocity.estimates <- function(emat,nmat,deltaT=1,smat=NULL,steady.state.cells=colnames(emat),kCells=10,cellKNN=NULL,kGenes=1,old.fit=NULL,mult=1e3,min.nmat.smat.correlation=0.05,min.nmat.emat.correlation=0.05, min.nmat.emat.slope=0.05, zero.offset=FALSE,deltaT2=1, fit.quantile=NULL, diagonal.quantiles=FALSE, show.gene=NULL, do.par=TRUE, cell.dist=NULL, emat.size=NULL, nmat.size=NULL, cell.emb=NULL, cell.colors=NULL, expression.gradient=NULL,residual.gradient=NULL, n.cores=defaultNCores(), verbose=TRUE) {
  if(!all(colnames(emat)==colnames(nmat))) stop("emat and nmat must have the same columns (cells)");
  if(!is.null(smat)) { if(!all(colnames(emat)==colnames(smat))) stop("smat must have the same columns (cells) as emat") }
  resl <- list();
  # bring matrices to the same gene set (just in case)
  vg <- intersect(rownames(emat),rownames(nmat));
  if(is.null(smat)) {
    emat <- emat[vg,]; nmat <- nmat[vg,]
  } else {
    vg <- intersect(vg,rownames(smat))
    emat <- emat[vg,]; nmat <- nmat[vg,]; smat <- smat[vg,]
  }
  if(!is.null(show.gene)) {
    if(!show.gene %in% rownames(emat)) { stop(paste("gene",show.gene,"is not present in the filtered expression matrices")) }
  }
  # TODO: add gene filtering options
  pcount <- 1;

  if(!is.null(cell.dist)) {
    if(class(cell.dist)!='dist') { stop("cell.dist must be of a class dist") }
    if(!all(labels(cell.dist)==colnames(emat))) {
      cat("matching cells between cell.dist and emat/nmat ... ")
      cell.dist <- as.matrix(cell.dist)
      cn <- intersect(colnames(emat),colnames(cell.dist))
      cell.dist <- as.dist(cell.dist[cn,cn]);
      emat <- emat[,cn]; nmat <- nmat[,cn];
      if(!is.null(smat)) { smat <- smat[,cn] }
      cat("done\n")
    }
  }
  
  # size estimates
  if(is.null(emat.size)) { emat.size <- Matrix::colSums(emat); }
  if(is.null(nmat.size)) { nmat.size <- Matrix::colSums(nmat); }
  emat.cs <- emat.size[colnames(emat)]/mult;
  nmat.cs <- nmat.size[colnames(nmat)]/mult;
  
  
  emat.log.norm <- log(as.matrix(t(t(emat)/emat.cs))+pcount);
  if(!is.null(old.fit)) { cellKNN <- old.fit[['cellKNN']]}
  knn.maxl <- 1e2
  if(kCells>1) {
    if(is.null(cellKNN)) {
      cat("calculating cell knn ... ")
      if(is.null(cell.dist)) {
        cellKNN <- balancedKNN(emat.log.norm,kCells,kCells*knn.maxl,n.threads=n.cores);
      } else {
        cellKNN <- balancedKNN(emat.log.norm,kCells,kCells*knn.maxl,n.threads=n.cores,dist=cell.dist);
      }
      diag(cellKNN) <- 1;
      resl$cellKNN <- cellKNN;
      cat("done\n")
    }
    rm(emat.log.norm);
    # smoothed matrices
    cat("calculating convolved matrices ... ")
    conv.emat <- emat %*% cellKNN[colnames(emat),colnames(emat)]
    conv.nmat <- nmat %*% cellKNN[colnames(nmat),colnames(nmat)]
    conv.emat.cs <- (emat.cs %*% cellKNN[colnames(emat),colnames(emat)])[1,]
    conv.nmat.cs <- (nmat.cs %*% cellKNN[colnames(nmat),colnames(nmat)])[1,]
    cat("done\n")
  } else {
    conv.emat <- emat; conv.nmat <- nmat; cellKNN <- NULL;
    conv.emat.cs <- emat.cs; conv.nmat.cs <- nmat.cs;
  }
  
  #browser()
  
  # size-normalized counts
  conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
  conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)

  # size-normalized counts
  emat.norm <- t(t(emat)/emat.cs)
  nmat.norm <- t(t(nmat)/nmat.cs)

  if(kGenes>1) {
    if(!is.null(old.fit) && !is.null(old.fit$geneKNN)) {
      geneKNN <- old.fit$geneKNN; 
    } else {
      cat("gene kNN ... ")
      geneKNN <- balancedKNN(t(log(as.matrix(conv.emat.norm)+pcount)),kGenes,kGenes*1.2e3,n.threads=n.cores); diag(geneKNN) <- 1;
    }
    resl$geneKNN <- geneKNN;


    # normalize contribution of different neighbor genes to match the median totals (to avoid distortions due to high-yielding genes)
    cat("scaling gene weights ... ")
    gt <- rowSums(conv.emat.norm)
    scaledGeneKNN <- t(apply(geneKNN,2,function(ii) pmin(1,median(gt[which(ii>0)])/gt) * ii))
    cat("convolving matrices ... ")
    conv.emat.norm <- scaledGeneKNN %*% conv.emat.norm;
    conv.nmat.norm <- scaledGeneKNN %*% conv.nmat.norm;
    
    cat("done\n")
  }
  
  if(!is.null(smat)) {
    
    if(kCells>1) {
      conv.smat <- smat %*% cellKNN[colnames(smat),colnames(smat)]
    } else {
      conv.smat <- smat
    }
    conv.smat.cs <- Matrix::colSums(conv.smat)/mult;
    conv.smat.norm <- t(t(conv.smat)/conv.smat.cs)

    if(kGenes>1) {
      conv.smat.norm <- scaledGeneKNN %*% conv.smat.norm;
    } 
    
    # use spanning reads to fit offset for the intronic reads, test correlation
    if(is.null(old.fit)) {
      cat("fitting smat-based offsets ... ")
      sfit <- data.frame(do.call(rbind,parallel::mclapply(sn(rownames(conv.emat.norm)),function(gn) {
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),s=conv.smat.norm[gn,steady.state.cells])
        sd <- lm(n~s,data=df)
        r <- with(df[df$s>0,],cor(n,s,method='spearman'),3)
        return(c(o=pmax(0,as.numeric(sd$coef[1])),s=as.numeric(sd$coef[2]),r=r))
      },mc.cores=n.cores,mc.preschedule=T)))
      cat("done\n")
      
    } else {
      sfit <- old.fit$sfit;
    }
  }

  resl$conv.nmat.norm <- conv.nmat.norm;
  resl$conv.emat.norm <- conv.emat.norm;

  # fit gamma, using the offset above
  if(!is.null(show.gene)) {
    gn <- show.gene;
    if(!is.null(cell.emb)) {
      # show embedding heatmaps
      cc <- intersect(rownames(cell.emb),colnames(conv.emat.norm));
      if(do.par) { par(mfrow=c(1,4), mar = c(2.5,2.5,2.5,0.5), mgp = c(1.5,0.65,0), cex = 0.85); }
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(conv.emat.norm[gn,cc],gradientPalette=expression.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'s'),axes=F); box();
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(conv.nmat.norm[gn,cc],gradientPalette=expression.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'u'),axes=F); box();
    }
    do <- NULL;
    if(!is.null(smat)) { # use smat-based offsets
      df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),o=sfit[gn,'o'])
      if(zero.offset) df$o <- 0;
    }  else { # calculate offset based on the nascent counts obsered for near-0 exonic levels
      df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]))
      o <- 0;
      df$o <- o;
      #zi <- emat[gn,steady.state.cells]==0;
      if(!zero.offset) { zi <- df$e<1/conv.emat.cs[steady.state.cells]; if(any(zi)) { o <- sum(df$n[zi])/(sum(zi)+1)} }
      df$o <- o;
      
      #table(zi)
      # if(any(zi)) {
      #   do <- lm(n~e,data=df[zi,])
      #   #summary(do)
      #   df$o <- max(0,do$coefficients[1])
      # }

      
    }
    
    
    #browser()
    d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
    cell.col <- ac(rep(1,nrow(df)),alpha=0.1); names(cell.col) <- rownames(df)
    if(!is.null(cell.colors)) { 
      cc <- intersect(names(cell.colors),rownames(df)); 
      cell.col[cc] <- cell.colors[cc]
    }
    plot(df$e,df$n,pch=21,bg=ac(cell.col,alpha=0.3),col=ac(1,alpha=0.1),cex=0.8,xlab='s',ylab='u',main=paste(gn,'fit'))
    if(!is.null(do)) {
      abline(do,lty=2,col=8)
    }
    
    # min/max fit
    if(!is.null(fit.quantile)) {
      if(diagonal.quantiles) {
        # determine maximum ranges 
        emax <- quantile(df$e,p=0.99)
        nmax <- quantile(df$n,p=0.99)
        if(emax==0) emax <- max(max(df$e),1e-3)
        if(nmax==0) nmax <- max(max(df$n),1e-3)
        x <- df$e/emax + df$n/nmax;
        eq <- quantile(x,p=c(fit.quantile,1-fit.quantile))
        if(!is.null(smat)) { # will use smat offset, so disregard lower quantile
          pw <- as.numeric(x>=eq[2])
        } else {
          pw <- as.numeric(x>=eq[2] | x<=eq[1])
        }
      } else {
        eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
        if(!is.null(smat) || zero.offset) { # will use smat offset, so disregard lower quantile
          pw <- as.numeric(df$e>=eq[2])
        } else {
          pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
        }
      }

      if(!is.null(smat) || zero.offset) { # use smat offset
        d <- lm(n~e+offset(o)+0,data=df,weights=pw);
      } else {
        d <- lm(n~e,data=df,weights=pw)
      }

      ## eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
      ## pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
      ## if(!is.null(smat)) { # use smat offset
      ##   d <- lm(n~e+offset(o)+0,data=df,weights=pw);
      ## } else {
      ##   d <- lm(n~e,data=df,weights=pw)
      ## }
    } 
    
    
    df <- df[order(df$e,decreasing=T),]; 
    lines(df$e,predict(d,newdata=df),lty=2,col=2)


    if(!is.null(cell.emb)) {
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(resid(d)[cc],gradientPalette=residual.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'resid'),axes=F); box();
    }
    if(kGenes>1) { return(invisible(geneKNN)) } else { return(1) }
  }
  
  cat("fitting gamma coefficients ... ")
  if(is.null(old.fit)) {
    ko <- data.frame(do.call(rbind,parallel::mclapply(sn(rownames(conv.emat.norm)),function(gn) {
      if(!is.null(smat)) { # use smat-based offsets
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),o=sfit[gn,'o'])
        if(zero.offset) df$o <- 0;
      }  else { # calculate offset based on the nascent counts obsered for near-0 exonic levels
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]))
        o <- 0;
        if(!zero.offset) { zi <- df$e<1/conv.emat.cs[steady.state.cells]; if(any(zi)) { o <- sum(df$n[zi])/(sum(zi)+1)} }
        df$o <- o;
      }
      if(is.null(fit.quantile)) {
        #d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
        d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
        return(c(o=df$o[1],g=as.numeric(coef(d)[1]),r=cor(df$e,df$n,method='spearman')))
      } else {
        if(diagonal.quantiles) {
          # determine maximum ranges 
          emax <- quantile(df$e,p=0.99)
          nmax <- quantile(df$n,p=0.99)
          if(emax==0) emax <- max(max(df$e),1e-3)
          if(nmax==0) nmax <- max(max(df$n),1e-3)
          x <- df$e/emax + df$n/nmax;
          eq <- quantile(x,p=c(fit.quantile,1-fit.quantile))
          if(!is.null(smat)) { # will use smat offset, so disregard lower quantile
            pw <- as.numeric(x>=eq[2])
          } else {
            pw <- as.numeric(x>=eq[2] | x<=eq[1])
          }
        } else {
          eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
          if(!is.null(smat) || zero.offset) { # will use smat offset, so disregard lower quantile
            pw <- as.numeric(df$e>=eq[2])
          } else {
            pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
          }
        }
        if(!is.null(smat) || zero.offset) { # use smat offset
          d <- lm(n~e+offset(o)+0,data=df,weights=pw);
          return(c(o=df$o[1],g=as.numeric(coef(d)[1]),r=cor(df$e,df$n,method='spearman')))
        } else {
          d <- lm(n~e,data=df,weights=pw)
          # note: re-estimating offset here
          return(c(o=as.numeric(coef(d)[1]),g=as.numeric(coef(d)[2]),r=cor(df$e,df$n,method='spearman')))
        }
      }
        
    },mc.cores=n.cores,mc.preschedule=T)))
    ko <- na.omit(ko)
    cat("done. succesfful fit for",nrow(ko),"genes\n")
  } else { full.ko <- ko <- na.omit(old.fit$ko); }

  if(!is.null(smat)) {
    sfit <- na.omit(sfit)
    ko <- ko[rownames(ko) %in% rownames(sfit),]; # omit genes for which sfit didn't work
    vi <- sfit$r > min.nmat.smat.correlation
    ko <- ko[vi,]
    if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-smat correlation\n")
  }
  
  full.ko <- ko;
  vi <- ko$r>min.nmat.emat.correlation
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat correlation\n")
  ko <- ko[vi,]
  
  vi <- ko$g>min.nmat.emat.slope
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat slope\n")
  ko <- ko[vi,]

  gamma <- ko$g; offset <- ko$o; names(gamma) <- names(offset) <- rownames(ko);
  cat("calculating RNA velocity shift ... ")
  if(kGenes>1) { # gene-convolved estimation
    # estimate M value
    npred <- gamma*conv.emat.norm[names(gamma),] + ko$o;
    npred[npred<0] <- 0;
    mval <- log2(conv.nmat.norm[names(gamma),]+pcount) - log2(npred+pcount);
    resl$mval <- mval;
    
    #resl$conv.deltaE <- t.get.projected.delta(conv.emat.norm,conv.nmat.norm,gamma,offset=offset,delta=deltaT)
    #resl$conv.projected <- t.get.projected.cell2(conv.emat.norm,emat.size,as.matrix(resl$conv.deltaE),mult = mult,delta=deltaT2);
    #resl$conv.projected[resl$conv.projected<0] <- 0;
    
    # switch back to non-gene-kNN conv.* matrices
    conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
    conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
    # estimate gamma
    cat("re-estimating gamma of individual genes ... ")
    
    am <- conv.nmat.norm[rownames(mval),]-offset; am[am<0] <- 0;
    fm <- log2(am) - mval - log2(conv.emat.norm[rownames(mval),])
    wm <- is.finite(fm)
    fm[!is.finite(fm)] <- 0;
    gammaA <- 2^(rowSums(fm * wm)/rowSums(wm))
    gammaA <- gammaA[is.finite(gammaA)];
 
    
    gamma <- gammaA;
    cat("done\n")

    
    # can estimate deltaE from the mval
    cat("calculating RNA velocity shift ... ")  
    # estimate delta from M value
    deltaE <- t.get.projected.delta.from.log2ratio(em=conv.emat.norm,gamma=gamma,r=mval,delta=deltaT)
    #deltaE <- t.get.projected.delta2(conv.emat.norm,conv.nmat,conv.nmat.cs,gamma,offset=offset,delta=deltaT)
  } else { # regular estimation
    #deltaE <- t.get.projected.delta2(conv.emat.norm,conv.nmat,conv.nmat.cs,gamma,offset=offset,delta=deltaT)
    deltaE <- t.get.projected.delta(conv.emat.norm,conv.nmat.norm,gamma,offset=offset,delta=deltaT)
  }

  resl$gamma <- gamma;
    
  cat("done\n")
  cat("calculating extrapolated cell state ... ")

 
  # reduced cell normalization (only genes for which momentum was estimated)
  emat.norm <- emat[rownames(emat) %in% rownames(deltaE),]
  #emat.sz <- Matrix::colSums(emat.norm)/mult;
  #browser()
  emat.sz <- emat.cs;
  emat.norm <- t(t(emat.norm)/(emat.sz));
  
  emn <- t.get.projected.cell2(emat.norm,emat.sz,as.matrix(deltaE),mult = mult,delta=deltaT2);
  #emn <- t.get.projected.cell(emat.norm,as.matrix(deltaE),target.mult = mult,model.mult=mult,delta=deltaT2,size.normalize=FALSE);
  
  #table(emat.norm[,cn]==0)
  #table(emn[,cn]==0)
  
  cat("done\n")
  full.ko$valid <- rownames(full.ko) %in% rownames(ko)
  resl <- c(resl,list(projected=emn,current=emat.norm,deltaE=deltaE,deltaT=deltaT,ko=full.ko,mult=mult,kCells=kCells));
  if(!is.null(smat)) { resl$sfit <- sfit }
  return(resl)
}

##' Structure-based gene velocity estimation
##'
##'
##' @param emat - spliced (exonic) count matrix
##' @param nmat - unspliced (nascent) count matrix
##' @param vel - initial gene-relative velocity estimates (output of the gene.relative.velocity.estimates function) 
##' @param base.df gene structure information data frame ($gene.df in output of read.gene.mapping.info()), containing the following columns ($il - total intronic length in log10(length+1) scale; $el - total exonic length; $nex - number of expressed (above some low threshold) exons; as well as optional $nipconc/$nipdisc giving number of concordant and discordant internal priming sites)
##' @param deltaT - amount of time to project the cell forward
##' @param smat - optional spanning read matrix (used in offset calculations)
##' @param kGenes - number of genes to use in evaluating trimmed mean of M values
##' @param kGenes.trim - number of genes to trim (from both ends)
##' @param smooth.kGenes - gene kNN pooling k value (used in the initial gene-relative fit)
##' @param kCells - number of k nearest neighbors (NN) to use in slope calculation smoothing
##' @param deltaT2 - scaling of the projected difference vector (normally should be set to 1)
##' @param min.gene.conuts - minimum number of spliced reads/molecules that a gene should have
##' @param min.gene.cells - minimum number of cells in which a gene should be expressed
##' @param min.intron.length - minimum exon length
##' @param min.exon.length - minimum exon length
##' @param top.global.pearson.deviance - maximum deviance threshold to filter out genes with very high unsplied counts (likely due to other processes)
##' @param cellKNN - optional pre-calculated cell KNN matrix
##' @param cell.dist - cell distance to use in cell kNN pooling calculations
##' @param fit.quantile perform gamma fit on a top/bottom quantiles of expression magnitudes
##' @param zero.offset force gene offsets to be zero (default if smat is not supplied), otherwise estimated from the lower quantile or quantile fit
##' @param diagonal.quantiles whether diagonal quantile determination should be used (if fit.quantile is specified)
##' @param m.pcount - pseudocount to be used in M value calculations (defaults to 5)
##' @param plot.model.fit plot gamma values predicted by the structure-bsaed model as a function of gene-relative gamma estimates.
##' @param n.cores - number of cores to use
##' @return a list with velocity results, including the current normalized expression state ($current), projected ($projected), unscaled transcriptional change ($deltaE), fit results ($ko, $sfit), optional cell pooling parameters ($cellKNN, $kCells), kNN-convolved normalized matrices (conv.nmat.norm and conv.emat.norm)
##' @examples
##' \dontrun{
##'  # emat / nmat are the spliced/unpsliced matrices respectively
##'  # rvel is a gene-relative velocity estimate
##'  # base.df (here dat$base.df) is a gene information table.
##'  #   For SMART-seq2, it is part of the \code{\link{read.smartseq2.bams}} output.
##'  #   For droplet data, this info can be obtained \code{\link{}}
##'  gvel <- global.velcoity.estimates(emat, nmat, rvel, dat$base.df, deltaT=1, kCells=5,
##'        kGenes = 15, kGenes.trim = 5, min.gene.cells = 0, min.gene.conuts = 500)
##'  
##' }
##' @export
global.velcoity.estimates <- function(emat,nmat,vel,base.df,deltaT=1,smat=NULL,kGenes=15,kGenes.trim=5,smooth.kGenes=0,kCells=10,deltaT2=1,min.gene.conuts=100,min.gene.cells=20,min.intron.length=10^3.5,min.exon.length=10^2.7,top.global.pearson.deviance=3,cellKNN=NULL,cell.dist=NULL,fit.quantile=NULL, zero.offset=NULL, diagonal.quantiles=FALSE, m.pcount=5,plot.model.fit=FALSE, n.cores=defaultNCores()) {

  if(is.null(zero.offset)) zero.offset <- is.null(smat); # set zero offset to true unless we have smat data
  
  mult <- vel$mult; # use the same library scale as in the supplied relative velocity estimates
  # reconsile gene lists
  gi <- intersect(intersect(rownames(base.df),rownames(emat)),rownames(nmat))
  emat <- emat[gi,]; nmat <- nmat[gi,]; 
  base.df <- base.df[gi,]

  # do some gene filtering
  vi <- rowSums(emat[rownames(base.df),])>min.gene.conuts & rowSums(emat[rownames(base.df),]>0)>min.gene.cells
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low emat levels\n")
  base.df <- base.df[vi,]
  
  vi <- base.df$il>log10(min.intron.length) & base.df$el>log10(min.exon.length) #requiring some minimum intronic and exonic lengths
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to insufficient exonic or intronic lengths\n")
  base.df <- base.df[vi,]

  mult <- vel$mult;
  
  # do a quick global model of total nascent reads as a function of total exonic reads and gene structural parameters to filter 
  # out the genes with very high nascent/exonic ratio - those are likely driven by other transcripts
  df <- data.frame(e=rowSums(emat[rownames(base.df),]), n=rowSums(nmat[rownames(base.df),]), base.df)
  #df$ir <- expr.lstat[rownames(df),'t']/expr.lstat[rownames(df),'i']
  df$eir <- base.df$il/base.df$el
  df$e <- log(df$e)
  gm <- MASS::glm.nb(n~.,data=df,link=log,init.theta=2)
  vi <- resid(gm,type='pearson') <= top.global.pearson.deviance;
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to excessive nascent counts\n")
  base.df <- base.df[vi,] 
  
  # reconcile gene lists and matrices
  gi <- intersect(rownames(base.df),rownames(emat))
  emat <- emat[gi,]; nmat <- nmat[gi,]; 
  base.df <- base.df[gi,]
  
  # start with gene-relative slopes
  gamma <- vel$gamma;
  gamma <- gamma[intersect(rownames(base.df),names(gamma))]
  #gamma <- ko$g; offset <- ko$o; names(gamma) <- names(offset) <- rownames(ko);
  cat("using relative slopes for",length(gamma),"genes to fit structure-based model ... ")
  
  df <- data.frame(k=log(gamma),base.df[names(gamma),]); # note we're working with log gamma to get good resolution of low values
  df$eir <- log(((10^df$il)-1)/((10^df$el)-1)) # intronic/exonic ratio
  df$e <- log(rowSums(emat[rownames(df),])) # total gene expression
  # genome-wide model fit
  # fit genome-wide model for the slope, based on the gene structural parameters and the total expression (exonic) magnitude
  if('nipconc' %in% colnames(base.df)) {
    cat("with internal priming info ... ")
    df <- cbind(df,data.frame(log10(base.df[names(gamma),c("nipconc","nipdisc")]+1))) 
    km <- mgcv::gam(k~s(il,e)+s(eir)+s(nex)+s(nipconc)+s(nipdisc),data=df,weights=sqrt(rowSums(emat[rownames(df),])))
  } else {
    km <- mgcv::gam(k~s(il,e)+s(eir)+s(nex),data=df,weights=sqrt(rowSums(emat[rownames(df),])))
  }
  cat(paste0(round((1-km$deviance/km$null.deviance)*100,1),"% deviance explained.\n"))

  if(plot.model.fit) {
    plot(df$k,predict(km),xlab=expression(paste('log[ gene-relative ',gamma,']')),ylab=expression(paste('log[ structure-based ',gamma,']')),pch=19,cex=0.7,col=ac(1,alpha=0.2));
    abline(a=0,b=1,col=2,lty=2)
    legend(x='bottomright',bty='n',legend=c(paste0(round((1-km$deviance/km$null.deviance)*100,1),"% deviance explained")))
  }
  
  
  # generate predictions for all genes
  df <- base.df; # note we're working with log gamma to get good resolution of low values
  df$eir <- log(((10^df$il)-1)/((10^df$el)-1)) # intronic/exonic ratio
  df$e <- log(rowSums(emat[rownames(df),])) # total gene expression
  cat("predicting gamma ... ")
  sqGammaPred <- predict(km,newdata=df,type='response')
  cat("done\n")
  
  # re-estimate offsets (and cellKNN) using relative fit
  cat("refitting offsets ... ")
  vel2 <- gene.relative.velocity.estimates(emat,nmat,smat=smat,kCells=kCells,kGenes=1,cell.dist=cell.dist,fit.quantile=fit.quantile,zero.offset=zero.offset,diagonal.quantiles=diagonal.quantiles)
  ko <- vel2$ko;
  cat("re-estimated offsets for",nrow(ko),"out of",nrow(emat),"genes\n")
  emat <- emat[rownames(ko),]; nmat <- nmat[rownames(ko),]
  sqGammaPred <- sqGammaPred[rownames(ko)]
  if(kCells>1) {
    cellKNN <- vel2$cellKNN;
    cat("calculating convolved matrices ... ")
    conv.emat <- emat %*% cellKNN[colnames(emat),colnames(emat)]
    conv.nmat <- nmat %*% cellKNN[colnames(nmat),colnames(nmat)]
    cat("done\n")    
  } else {
    conv.emat <- emat; conv.nmat <- nmat;
  }
  
  # size estimates

  conv.emat.cs <- Matrix::colSums(conv.emat)/mult;
  conv.nmat.cs <- Matrix::colSums(conv.nmat)/mult;
  # size-normalized counts
  conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
  conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
  

  
  # TODO: we want to restrict the set of genes that can be neighbors
  cat("calculating gene knn ... ")
  emm <- log10(as.matrix(conv.emat.norm)+1)
  gknn <- balancedKNN(t(emm),kGenes,nrow(emm),n.threads=n.cores); diag(gknn) <- 1;
  cat("done\n")
  cat("estimating M values ... ")
  # amount of nascent transcription predicted by the model-based k esimates
  npred <- t(t(conv.emat.norm[names(sqGammaPred),]*exp(as.numeric(sqGammaPred)))*conv.nmat.cs)
  # adjust for offset
  npred <- npred+ko$o
  
  mval <- log2((conv.nmat[rownames(npred),]+m.pcount)/(npred+m.pcount))
  cat("adjusting mval offsets ... ")
  # estimate median M value (log2 observed nascent / expected nascent ratio) across kNN genes
  mmval <- do.call(rbind,parallel::mclapply(sn(colnames(gknn)),function(gn) {
    gin <- names(which(gknn[,gn]>0)) # neighbor gene names
    #x <- apply(mval[gin,],2,median) # works well too
    x <- apply(mval[gin,],2,mean,trim=kGenes.trim)
  },mc.cores=n.cores,mc.preschedule=TRUE))
  
  if(kGenes>1) { # gene-convolved estimation
    # switch back to non-gene-kNN conv.* matrices
    conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
    conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
  }

  # adjust gamma predictions
  cat("re-estimating gamma ... ")

  offset <- ko[,'o']

  ### alternative gamma fit procedure
  ## gammaA <- unlist(mclapply(sn(rownames(mval)),function(gn) {
  ##   # here we try to optimize k to match the expression change based on mmval
  ##   df <- data.frame(n=conv.nmat.norm[gn,],e=conv.emat.norm[gn,],m=2^mval[gn,],o=ko[gn,'o'])
  ##   df$na <- df$n-df$o; df$na[df$na<0] <- 0; # apply offset
  ##   pc <- mean(df$e)/5; # min count
  ##   # fitness function, capturing discrepancy in deltaE
  ##   #f <- function(k) { ep <- (df$e+pc)*df$m - df$e; np <- df$n/k-df$e; sum(sqrt(abs(ep-np))) }
  ##   f <- function(k) { egt <- exp(-k); ep <- (df$e+pc)*(egt*(1-df$m)+df$m) - df$e; np <- df$e*egt + (df$na)/k*(1-egt)-df$e; sum(sqrt(abs(ep-np))) }
  ##   # poor man's optimization here, to guide interval methods
  ##   iv <- 10^seq(-4,2,by=0.1)
  ##   ivv <- unlist(lapply(iv,f))
  ##   mi <- which.min(ivv)
  ##   ov <- optim(iv[mi],f,lower=iv[max(1,mi-2)],upper=iv[min(length(iv),mi+2)],method='L-BFGS-B')
  ##   if(ov$value<ivv[mi]) { return((ov$par))} else { return((iv[mi]))}
  ## },mc.cores=n.cores,mc.preschedule=T))


  am <- conv.nmat.norm[rownames(mmval),]-offset; am[am<0] <- 0;
  fm <- log2(am) - mmval - log2(conv.emat.norm[rownames(mmval),])
  wm <- is.finite(fm)
  fm[!is.finite(fm)] <- 0;
  gammaA <- 2^(rowSums(fm * wm)/rowSums(wm))
  gammaA <- gammaA[is.finite(gammaA)];

  #plot(log(old.gammaA),log(gammaA)); abline(a=0,b=1,lty=2,col=2);
  
  cat("done\n")
  
  # can estimate deltaE from the mval
  cat("calculating RNA velocity shift ... ")  
  deltaE <- t.get.projected.delta.from.log2ratio(em=conv.emat.norm,gamma=gammaA,r=mmval,delta=deltaT)
  cat("done\n")
  cat("calculating extrapolated cell state ... ")
  emat.size <- Matrix::colSums(emat)/mult;
  
  em <- as.matrix(t(t(emat)/emat.size))
  em <- em[rownames(em) %in% rownames(deltaE),]
  emn <- t.get.projected.cell2(em,emat.size,as.matrix(deltaE),mult = mult,delta=deltaT2);
  cat("done\n")
  resl <- list(projected=emn,current=em,deltaE=deltaE,deltaT=deltaT,mval=mval,mult=mult,vel2=vel2,gammaA=gammaA);
  return(resl)
  
}

##' Filter genes by requirining minimum average expression within at least one of the provided cell clusters
##'
##' @param emat spliced (exonic) count matrix
##' @param clusters named cell factor defining clusters
##' @param min.max.cluster.average required minimum average expression count (no normalization is perfomed)
##' @return filtered emat matrix
##' @export
filter.genes.by.cluster.expression <- function(emat,clusters,min.max.cluster.average=0.1) {
  if(!any(colnames(emat) %in% names(clusters))) stop("provided clusters do not cover any of the emat cells!")
  vc <- intersect(colnames(emat),names(clusters))
  cl.emax <- apply(do.call(cbind,tapply(vc,as.factor(clusters[vc]),function(ii) Matrix::rowMeans(emat[,ii]))),1,max)
  vi <- cl.emax>min.max.cluster.average;
  emat[vi,]
}

##' PCA-based visualization of the velocities
##'
##' @param vel velocity estimation (gene-relative or global)
##' @param nPcs number of successive PCs to visualize
##' @param cell.colors a named vector of cell colors for visualization
##' @param scale scale to use for expression state transform (default: 'log', other possible values are 'sqrt','linear')
##' @param plot.cols number of columns into which to arrange the plots
##' @param norm.nPcs optional total number of PCs to use for velocity magnitude normalization
##' @param do.par whether to set up graphical parameters of a plot
##' @param pc.multipliers an optional vector multipliers for the cell PC scores (useful for reorienting the PCs)
##' @param show.grid.flow whether a grid flow should be shown
##' @param grid.n number of grid points (on each axis)
##' @param grid.sd standard deviation of the grid
##' @param arrow.scale scale multiplier for the velocity estimates
##' @param min.grid.cell.mass minimum cellular mass
##' @param min.arrow.size minimum size of an arrow to show
##' @param pcount pseudocount
##' @param arrow.lwd thickness of arrows to plot
##' @param size.norm whether to rescale current and projected states by cell size (default=FALSE)
##' @param return.details whether to return detailed output
##' @param plot.grid.points whether to show dots at every grid point
##' @param fixed.arrow.length whether to use fixed-size arrow
##' @param max.grid.arrow.length limit to the size of the arrows that could be shown (when fixed.arrow.length=FALSE)
##' @param n.cores number of cores to use in the calculations
##' @param ... extra parameters are passed to plot() function
##' @return If return.details=F, returns invisible list containing PCA info (epc) and projection of velocities onto the PCs (delta.pcs). If return.details=T, returns an extended list that can be passed into p1 app for velocity visualization.
##' @export
pca.velocity.plot <- function(vel,nPcs=4,cell.colors=NULL,scale='log',plot.cols=min(3,nPcs-1),norm.nPcs=NA,do.par=T, pc.multipliers=NULL, show.grid.flow=FALSE, grid.n=20, grid.sd=NULL, arrow.scale=1, min.grid.cell.mass=1, min.arrow.size=NULL, pcount=1, arrow.lwd=1, size.norm=FALSE, return.details=FALSE, plot.grid.points=FALSE, fixed.arrow.length=FALSE,max.grid.arrow.length=NULL, n.cores=defaultNCores(), ...) {
  x0 <- vel$current;
  x1 <- vel$projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(x0)),alpha=0.3); names(cell.colors) <- colnames(x0) }
  # rescale to the same size
  if(size.norm) {
    cat("rescaling ... ")
    sz <- Matrix::colSums(x0);
    x0 <- t(t(x0)/sz)*mean(sz)
    x1 <- t(t(x1)/Matrix::colSums(x1))*mean(sz);
  }
  # transform
  if(scale=='log') { 
    cat("log ... ")
    x0.log <- log2(x0+pcount)
    x1.log <- log2(x1+pcount)
  } else if(scale=='sqrt') {
    cat("sqrt ... ")
    x0.log <- sqrt(x0)
    x1.log <- sqrt(x1)
  } else { # linear
    cat("linear ... ")
    x0.log <- x0
    x1.log <- x1
  }
  
  cat("pca ... ")
  cent <- rowMeans(x0.log);
  epc <- pcaMethods::pca(t(x0.log-cent),center=F,nPcs=ifelse(is.na(norm.nPcs),nPcs,norm.nPcs))
  
  if(!is.null(pc.multipliers)) { # apply multipliers (used for flipping the direction of PCs in the plots)
    if(length(pc.multipliers)!=nPcs) stop("pc.multipliers must be a vector equal in length to the number of PCs")
    cat("pc multipliers ... ")
    epc@loadings <- t(t(epc@loadings)*pc.multipliers)
    epc@scores <- scale(epc@completeObs,scale=F,center=T) %*% epc@loadings;
  }
  
  x1.scores <- t(x1.log - cent) %*% epc@loadings
  
  # normalize velocities ...?
  cat("delta norm ... ")
  delta.pcs <- as.matrix(x1.scores-epc@scores)
  if(!is.na(norm.nPcs)) {
    delta.pcs <- delta.pcs/mean(sqrt(rowSums(delta.pcs^2))) # suggested by Gioele, unsure about this ....
  }
  
  # browser()
  # z <- as.matrix(t(x1.log-x0.log)) %*% epc@loadings
  # 
  # hist(apply(vel$deltaE,2,mean))
  # summary(apply(vel$deltaE,2,mean))
  # hist(apply(as.matrix(x1.log-x0.log),2,mean))
  # summary(apply(as.matrix(x1.log-x0.log),2,mean))
  # 
  # z <- t(as.matrix(vel$deltaE)[rownames(epc@loadings),]) %*%  epc@loadings
  # str(z)
  # str(delta.pcs)
  # cn <- 'L6'
  # cn <- 'O15'
  # z <- (x1.log-x0.log)[,cn] * epc@loadings[,3]
  # z2 <- vel$deltaE[names(z),cn] * epc@loadings[,3]
  # sort(z,d=T)[1:20]
  # sum(z)
  # delta.pcs[cn,3]
  # summary(delta.pcs[,3])
  # str(epc@loadings[,3])
  # sort(delta.pcs[,3],d=T)[1:10]
  # z <- rowMeans(x1.log-x0.log) * epc@loadings[,3]
  #summary(z)
  #sort(z,d=T)[1:10]

  delta.pcs <- delta.pcs *arrow.scale;
  cat("done\n")
  if(do.par) par(mfrow=c(ceiling((nPcs-1)/plot.cols),plot.cols), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  vinfo <- lapply(1:(nPcs-1),function(i) {
    pos <- epc@scores[,c((i-1)+1,(i-1)+2)];
    #ppos <- x1.scores[,c((i-1)+1,(i-1)+2)];
    ppos <- pos+delta.pcs[,c((i-1)+1,(i-1)+2)];
    plot(pos,bg=cell.colors[rownames(pos)],pch=21,col=ac(1,alpha=0.3),lwd=0.5,xlab=paste("PC",(i-1)+1),ylab=paste("PC",(i-1)+2),axes=T,main=paste('PC',(i-1)+1,' vs. PC',(i-1)+2,sep=''),  ...); box();
    
    if(show.grid.flow) { # show grid summary of the arrows
      # arrow estimates for each cell
      ars <- data.frame(pos[,1],pos[,2],ppos[,1],ppos[,2])
      colnames(ars) <- c('x0','y0','x1','y1')
      arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
      rownames(ars) <- rownames(arsd) <- rownames(pos);
      
      # set up a grid
      rx <- range(c(range(ars$x0),range(ars$x1)))
      ry <- range(c(range(ars$y0),range(ars$y1)))
      gx <- seq(rx[1],rx[2],length.out=grid.n)
      gy <- seq(ry[1],ry[2],length.out=grid.n)
      
      # for each grid point calculate Gaussian-weighted delta average
      if(is.null(grid.sd)) {
        grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
        cat("grid.sd=",grid.sd," ")
      }

      if(is.null(min.arrow.size)) {
        min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
        cat("min.arrow.size=",min.arrow.size," ")
      }

      if(is.null(max.grid.arrow.length)) {
        max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
        cat("max.grid.arrow.length=",max.grid.arrow.length," ")
      }


      garrows <- do.call(rbind,lapply(gx,function(x) {
        # cell distances (rows:cells, columns: grid points)
        cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
        cw <- dnorm(cd,sd=grid.sd)
        # calculate x and y delta expectations
        gw <- Matrix::colSums(cw)
        cws <- pmax(1,Matrix::colSums(cw));
        gxd <- Matrix::colSums(cw*arsd$xd)/cws
        gyd <- Matrix::colSums(cw*arsd$yd)/cws
        
        al <- sqrt(gxd^2+gyd^2);
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
        
        cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])
      }))
      colnames(garrows) <- c('x0','y0','x1','y1')
      
      # plot
      if(fixed.arrow.length) {
        suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=0.05,lwd=arrow.lwd))
      } else {
        alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
        # can't specify different arrow lengths in one shot :(
        suppressWarnings(lapply(1:nrow(garrows),function(i) arrows(garrows[i,1],garrows[i,2],garrows[i,3],garrows[i,4],length=alen[i],lwd=arrow.lwd)))
      }
      if(plot.grid.points) points(rep(gx,each=length(gy)),rep(gy,length(gx)),pch='.',cex=1e-1,col=ac(1,alpha=0.4))
    
      if(return.details) { # for the p1 app
        # calculate expression shift
        cat("expression shifts .")
        # for individual cells
        es <- as.matrix(epc@loadings[,c((i-1)+1,(i-1)+2)] %*% t(delta.pcs[,c((i-1)+1,(i-1)+2)]))
        
        cat(".");
        gs <- epc@loadings[,c((i-1)+1,(i-1)+2)] %*% rbind(garrows[,3]-garrows[,1],garrows[,4]-garrows[,2])

        # note: here we're using deltaE vector, which may be normalized a bit differently from the $current/$projectted that was used above
        nd <- as.matrix(vel$deltaE)
        if(scale=='log') {
          nd <- (log10(abs(nd)+1)*sign(nd))
        } else if(scale=='sqrt') {
          nd <- (sqrt(abs(nd))*sign(nd))
        }
        cat(".");
        # velocity for the grid (weight-averaged velocity vectors)
        
        
        gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
          cw <- dnorm(cd,sd=grid.sd)
          # calculate x and y delta expectations
          gw <- Matrix::colSums(cw)
          cws <- pmax(1,Matrix::colSums(cw));
          cw <- t(t(cw)/cws)
          gxd <- Matrix::colSums(cw*arsd$xd)
          gyd <- Matrix::colSums(cw*arsd$yd)
          al <- sqrt(gxd^2+gyd^2);
          vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
          if(any(vg)) {
            z <- nd %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        cat(". done\n")
        
        return(invisible(list(garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,gvel=gv,geshifts=gs,scale=scale,emb=pos,epc=epc)))
      }
      
    } else {
      # draw individual arrows
      grid();
      suppressWarnings(arrows(pos[,1],pos[,2],ppos[,1],ppos[,2],length=0.05,lwd=arrow.lwd))
    }
  })
  cat("done\n")
  if(return.details) { return(vinfo) }
  return(invisible(list(epc=epc,delta.pcs=delta.pcs)))
  
}
##' Joint t-SNE visualization of the velocities by joint t-SNE embedding of both current and extraploated cell positions
##'
##' @param vel velocity result
##' @param cell.colors named color vector for the cells
##' @param scale whether to rescale current/projected
##' @param do.par whether to reset par (default=T)
##' @param delta.norm whether to renormalize velocities following PCA projection
##' @param nPcs number of PCs onto which project the velocities
##' @param norm.nPcs number of PCs to use for velocity normalization
##' @param perplexity perplexity parameter to use in joint t-SNE calculation
##' @param show.grid.flow whether grid flow pattern should be drawn
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation of each grid point (used to determine the averaging radius for each grid point)
##' @param min.grid.cell.mass minimal number of cells around a grid point required for the grid point to show up
##' @param pcount pseudocount
##' @param verbose whether to show messages
##' @param min.arrow.median.ratio minimal ratio of arrow length (to the median arrow length) below which the arrows are not drawn (default=1/10)
##' @param max.arrow.quantile max arrow quantile that's used for arrow size calculation (default=0.9)
##' @param arrow.scale scaling factor for the arrows
##' @param arrow.lwd arrow line width
##' @param xlab x axis label
##' @param ylab y axis label
##' @param n.cores number of cores to use
##' @param size.norm whether to re-normalize current and projected cell sizes
##' @param ... extra parameters are passed to plot() routine.
##' @return invisible list containing embedding positions of current state (current.emb) and extrapolated states (projected.emb)
##' @export
tSNE.velocity.plot <- function(vel,cell.colors=NULL,scale='log',do.par=T, delta.norm=TRUE, nPcs=15, norm.nPcs=nPcs*10, perplexity=ncol(vel$current)/3, show.grid.flow=FALSE, grid.n=20, grid.sd=NULL, min.grid.cell.mass=1, pcount=0.1, verbose=TRUE, min.arrow.median.ratio=1/10, max.arrow.quantile=0.9, arrow.scale=1, arrow.lwd=1, xlab="", ylab="", n.cores=defaultNCores(), size.norm=TRUE, ...) {
  x0 <- vel$current;
  x1 <- vel$projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(x0)),alpha=0.3); names(cell.colors) <- colnames(x0) }
  
  if(size.norm) {
    # rescale to the same size
    cat("rescaling ... ")
    sz <- Matrix::colSums(x0);
    x0 <- t(t(x0)/sz)*mean(sz)
    x1 <- t(t(x1)/Matrix::colSums(x1))*mean(sz);
  }
  # transform
  if(scale=='log') { 
    cat("log ... ")
    x0.log <- log2(x0+pcount)
    x1.log <- log2(x1+pcount)
  } else if(scale=='sqrt') {
    cat("sqrt ... ")
    x0.log <- sqrt(x0)
    x1.log <- sqrt(x1)
  } else { # linear
    cat("linear ... ")
    x0.log <- x0
    x1.log <- x1
  }
  if(!is.null(nPcs)) { # reduce using PCA first
    cat("pca ... ")
    cent <- rowMeans(x0.log);
    epc <- pcaMethods::pca(t(x0.log-cent),center=F,nPcs=ifelse(is.na(norm.nPcs),nPcs,norm.nPcs))
    x0.log <- epc@scores;
    x1.log <- t(x1.log - cent) %*% epc@loadings
    if(delta.norm) {
      # normalize velocities ...?
      cat("delta norm ... ")
      delta.pcs <- x1.log-x0.log;
      if(!is.na(norm.nPcs)){ 
        delta.pcs <- delta.pcs/mean(sqrt(rowSums(delta.pcs^2))) # ? unsure about this (cell-wise L2 norm)
      }
      x1.log <- x0.log+delta.pcs;
    }
    # drop extra Pcs 
    x0.log <- t(x0.log[,1:nPcs])
    x1.log <- t(x1.log[,1:nPcs])
  }
  cat("tSNE ...")
  emb <- Rtsne::Rtsne(t(cbind(as.matrix(x0.log),as.matrix(x1.log))), num_threads=n.cores, perplexity=perplexity, verbose=verbose)$Y;
  x0.emb <- emb[1:ncol(x0.log),]
  x1.emb <- emb[-(1:ncol(x0.log)),]
  rownames(x0.emb) <- rownames(x1.emb) <- colnames(x0.log);
  
  cat("delta norm ... ") # again, somewhat unsure about this part
  delta.emb <- x1.emb - x0.emb;
  asize <- rowSums(delta.emb^2);
  no.arrow <- asize<= median(asize)*min.arrow.median.ratio;
  # restrict size top top 90% quantile
  delta.emb <- delta.emb/asize * pmin(asize,2*quantile(asize,p=max.arrow.quantile))*arrow.scale;
  delta.emb[no.arrow,] <- 0;
  cat("done\n") 
  if(do.par) par(mfrow=c(1,1), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  plot(x0.emb,bg=cell.colors[rownames(x0.emb)],pch=21,col=ac(1,alpha=0.3),xlab=ylab,ylab=xlab, ... ); box();
  
  
  if(show.grid.flow) { # show grid summary of the arrows
    # arrow estimates for each cell
    ars <- data.frame(x0.emb[,1],x0.emb[,2],x0.emb[,1]+delta.emb[,1],x0.emb[,2]+delta.emb[,2])
    colnames(ars) <- c('x0','y0','x1','y1')
    arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
    
    # set up a grid
    cat("grid estimates ... ")
    rx <- range(c(range(ars$x0),range(ars$x1)))
    ry <- range(c(range(ars$y0),range(ars$y1)))
    gx <- seq(rx[1],rx[2],length.out=grid.n)
    gy <- seq(ry[1],ry[2],length.out=grid.n)
    
    # for each grid point calculate Gaussian-weighted delta average
    if(is.null(grid.sd)) {
      grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
    }
    ginfo <- lapply(gx,function(x) {
      # cell distances (rows-cells,columsn - grid points)
      cd <- sqrt(outer(x0.emb[,2],-gy,'+')^2 + (x-x0.emb[,1])^2)
      cw <- dnorm(cd,sd=grid.sd)
      # calculate x and y delta expectations
      gw <- Matrix::colSums(cw)
      cws <- pmax(1,Matrix::colSums(cw));
      gxd <- Matrix::colSums(cw*arsd$xd)/cws
      gyd <- Matrix::colSums(cw*arsd$yd)/cws
      
      vg <- gw>=min.grid.cell.mass
      if(any(vg)) {
        suppressWarnings(arrows(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg],length=0.05,lwd=arrow.lwd))
      }
      points(rep(x,length(gy)),gy,pch='.')
    })
    cat("done\n")
  } else {
    # draw individual arrows
    suppressWarnings(arrows(x0.emb[,1],x0.emb[,2],x0.emb[,1]+delta.emb[,1],x0.emb[,2]+delta.emb[,2],length=0.05,lwd=arrow.lwd))
  }

  return(invisible(list(current.emb=x0.emb,projected.emb=x1.emb+delta.emb)))
}

##' Visualize RNA velocities on an existing embedding using correlation-based transition probability matrix within the kNN graph
##'
##' @param emb embedding onto which to project the velocities; The dimensions of coordinates should be on the order of 10x10 for the default values to make sense.
##' @param vel velocity estimates (e.g. returned by gene.relative.velocity.estimates() )
##' @param n neighborhood size (default=100 cells)
##' @param cell.colors name vector of cell colors
##' @param corr.sigma sigma parameter used to translate velocity-(expression delta) correlation into a transition probability
##' @param show.grid.flow whether to show grid velocity summary
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation (in embedding coordinate space) used to determine the weighting of individual cells around each grid point
##' @param min.grid.cell.mass minimal cell "mass" (weighted number of cells) around each grid point required for it to show up
##' @param min.arrow.size minimal arrow size
##' @param arrow.scale arrow scale multiplier
##' @param max.grid.arrow.length minimal arrow size
##' @param fixed.arrow.length whether to use fixed arrow width (default=FALSE)
##' @param plot.grid.points whether to mark all grid points with dots (even if they don't have valid velocities)
##' @param scale velocity scale to use (default: 'log', other values: 'sqrt','rank','linear')
##' @param nPcs number of PCs to use for velocity regularization (default NA, turns off regularization)
##' @param arrow.lwd arrow width (under fixed.arrow.length=T)
##' @param xlab x axis label
##' @param ylab y axls label
##' @param n.cores number of cores to use
##' @param do.par whether to reset plotting parameters
##' @param show.cell whether to show detailed velocity estimates for a specified cell
##' @param cell.border.alpha transparency for the cell border
##' @param cc velocity-(exprssion delta) correlation matrix (can be passed back from previous results, as $cc) to save calculation time when replotting the same velocity estimates on the same embedding with different parameters
##' @param return.details whether to return detailed output (which can be passed to p1 app for visualization)
##' @param expression.scaling whether to scale the velocity length by the projection of velocity onto the expected expression change (based on the transition probability matrix)
##' @param ... extra parameters passed to plot() function
##' @return if return.details=F, returns invisible list containing transition probability matrix ($tp) and the velocity-(expression delta) correlation matrix ($cc). If return.details=T, returns a more extended list that can be passed as veloinfo to pagoda2::p2.make.pagoda1.app() for visualization
##' @export
show.velocity.on.embedding.cor <- function(emb,vel,n=100,cell.colors=NULL, corr.sigma=0.05, show.grid.flow=FALSE, grid.n=20, grid.sd=NULL, min.grid.cell.mass=1, min.arrow.size=NULL, arrow.scale=1, max.grid.arrow.length=NULL, fixed.arrow.length=FALSE, plot.grid.points=FALSE, scale='log', nPcs=NA,  arrow.lwd=1, xlab="", ylab="", n.cores=defaultNCores(), do.par=T, show.cell=NULL, cell.border.alpha=0.3,cc=NULL, return.details=FALSE, expression.scaling=FALSE,  ...) {
  randomize <- FALSE;
  if(do.par) par(mfrow=c(1,1), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  celcol <- 'white'
  if(is.null(show.cell)) { celcol <- cell.colors[rownames(emb)] }
  plot(emb,bg=celcol,pch=21,col=ac(1,alpha=cell.border.alpha), xlab=xlab, ylab=ylab, ...);
  
  #plot(emb,bg=cell.colors[rownames(emb)],pch=21,col=ac(1,alpha=0.3), xlab=xlab, ylab=ylab);
  em <- as.matrix(vel$current); 
  ccells <- intersect(rownames(emb),colnames(em));
  em <- em[,ccells]; emb <- emb[ccells,]
  nd <- as.matrix(vel$deltaE[,ccells])
  cgenes <- intersect(rownames(em),rownames(nd));
  nd <- nd[cgenes,]; em <- em[cgenes,]
  #browser()
  if(randomize) {
    # randomize cell and sign for each gene
    nd <- t(apply(nd,1,function(x) (rbinom(length(x),1,0.5)*2-1)*abs(sample(x))))
  }
  #vg <- rownames(em) %in% rownames(r)

  
  if(is.null(cc)) {
    # cosine projections
    cat("delta projections ... ")

    if(scale=='log') {
      cat("log ")
      cc <- colDeltaCorLog10(em,(log10(abs(nd)+1)*sign(nd)),nthreads=n.cores);
    } else if(scale=='sqrt') {
      cat("sqrt ")
      cc <- colDeltaCorSqrt(em,(sqrt(abs(nd))*sign(nd)),nthreads=n.cores);
    } else if(scale=='rank') {
      cat("rank ")
      cc <- colDeltaCor((apply(em,2,rank)),(apply(nd,2,rank)),nthreads=n.cores);
    } else { # linear
      cat("linear ")
      cc <- colDeltaCor(em,nd,nthreads=n.cores);
    }
    colnames(cc) <- rownames(cc) <- colnames(em)
    diag(cc) <- 0;
  }

  cat("knn ... ")
  if(n>nrow(cc)) { n <- nrow(cc) }
  # TODO: add kNN based on high-dimensional correlation or Euclidean distances
  # define kNNs based on the embedding (L2 distance)
  emb.knn <- balancedKNN(t(emb),k=n,maxl=nrow(emb),dist='euclidean',n.threads=n.cores)
  diag(emb.knn) <- 1
  # caluclate transition probabilities (from col to row)
  cat("transition probs ... ")
  tp <- exp(cc/corr.sigma)*emb.knn
  #diag(tp) <- 0; #  should we allow the self-corelation for scaling?
  tp <- t(t(tp)/Matrix::colSums(tp)); # tp shows transition from a given column cell to different row cells
  tp <- as(tp,'dgCMatrix')
  cat("done\n")
  if(!is.null(show.cell)) {
    i <- match(show.cell,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # plot transition prob for a given cell
    points(emb,pch=19,col=ac(val2col(tp[rownames(emb),show.cell],gradient.range.quantile=1),alpha=0.5))
    points(emb[show.cell,1],emb[show.cell,2],pch=3,cex=1,col=1)
    di <- t(t(emb)-emb[i,])
    di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
    dir <- Matrix::colSums(di*tp[,i]) 
    dic <- Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center
    dia <- dir-dic;
    #browser()
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],length=0.05,lwd=1,col='blue'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,col='red'))
    suppressWarnings(arrows(emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,lty=1,col='grey50'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dia[1],emb[colnames(em)[i],2]+dia[2],length=0.05,lwd=1,col='black'))
  } else {
    # arrow estimates for each cell
    cat("calculating arrows ... ")
    arsd <- data.frame(t(embArrows(emb,tp,arrow.scale,n.cores)))
    rownames(arsd) <- rownames(emb);
    
    if(expression.scaling) {
      tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
      es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
      # project velocity onto expression shift
      #pm <- as.matrix(t(vel$deltaE)/sqrt(colSums(vel$deltaE*vel$deltaE)))[colnames(es),] * (t(es)/sqrt(colSums(es*es)))
      #pl <- pmax(0,apply(pm,1,sum))
      pl <- pmin(1,pmax(0,apply(as.matrix(vel$deltaE[,colnames(es)]) * es, 2, sum)/sqrt(colSums(es*es))))
      
      
      arsd <- arsd * pl;
    }
    
    
    ars <- data.frame(cbind(emb,emb+arsd));
    colnames(ars) <- c('x0','y0','x1','y1')
    colnames(arsd) <- c('xd','yd')
    rownames(ars) <- rownames(emb);
    cat("done\n")
    
    
    if(show.grid.flow) { # show grid summary of the arrows

    # set up a grid
    cat("grid estimates ... ")
    rx <- range(c(range(ars$x0),range(ars$x1)))
    ry <- range(c(range(ars$y0),range(ars$y1)))
    gx <- seq(rx[1],rx[2],length.out=grid.n)
    gy <- seq(ry[1],ry[2],length.out=grid.n)
    
    # for each grid point calculate Gaussian-weighted delta average
    if(is.null(grid.sd)) {
      grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
      cat("grid.sd=",grid.sd," ")
    }
    if(is.null(min.arrow.size)) {
      min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
      cat("min.arrow.size=",min.arrow.size," ")
    }
    if(is.null(max.grid.arrow.length)) {
      max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
      cat("max.grid.arrow.length=",max.grid.arrow.length," ")
    }
    
    garrows <- do.call(rbind,lapply(gx,function(x) {
      # cell distances (rows:cells, columns: grid points)
      cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
      cw <- dnorm(cd,sd=grid.sd)
      # calculate x and y delta expectations
      gw <- Matrix::colSums(cw)
      cws <- pmax(1,Matrix::colSums(cw));
      gxd <- Matrix::colSums(cw*arsd$xd)/cws
      gyd <- Matrix::colSums(cw*arsd$yd)/cws
      
      al <- sqrt(gxd^2+gyd^2);
      vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
      
      cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])
    }))
    colnames(garrows) <- c('x0','y0','x1','y1')

    # plot
    if(fixed.arrow.length) {
      suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=0.05,lwd=arrow.lwd))
    } else {
      alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
      # can't specify different arrow lengths in one shot :(
      #suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=alen,lwd=arrow.lwd))
      suppressWarnings(lapply(1:nrow(garrows),function(i) arrows(garrows[i,1],garrows[i,2],garrows[i,3],garrows[i,4],length=alen[i],lwd=arrow.lwd)))
    }
    if(plot.grid.points) points(rep(gx,each=length(gy)),rep(gy,length(gx)),pch='.',cex=1e-1,col=ac(1,alpha=0.4))
    
    cat("done\n")
    
    if(return.details) { # for the p1 app
      # calculate expression shift
      cat("expression shifts .")
      # for individual cells
      
      scale.int <- switch(scale,'log'=2,'sqrt'=3,1)
      #es <- expectedExpressionShift(e=as.matrix(em),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
      if(!expression.scaling) { #otherwise it has already been calculated
        tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
        #es <- expectedExpressionShift(e=as.matrix(em %*% as.matrix(tpb)),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
        es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
      }
      cat(".");
      # for the grid
      gs <- do.call(cbind,parallel::mclapply(gx,function(x) {
        # cell distances (rows:cells, columns: grid points)
        cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
        cw <- dnorm(cd,sd=grid.sd)
        # calculate x and y delta expectations
        gw <- Matrix::colSums(cw)
        cws <- pmax(1,Matrix::colSums(cw));
        cw <- t(t(cw)/cws)
        gxd <- Matrix::colSums(cw*arsd$xd)
        gyd <- Matrix::colSums(cw*arsd$yd)
        al <- sqrt(gxd^2+gyd^2);
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
        if(any(vg)) {
          z <- es %*% cw[,vg]
        } else { NULL }
      },mc.cores=n.cores,mc.preschedule=T))

      if(scale=='log') {
        nd <- (log10(abs(nd)+1)*sign(nd))
      } else if(scale=='sqrt') {
        nd <- (sqrt(abs(nd))*sign(nd))
      }
      cat(".");
      # velocity for the grid
      gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
        # cell distances (rows:cells, columns: grid points)
        cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
        cw <- dnorm(cd,sd=grid.sd)
        # calculate x and y delta expectations
        gw <- Matrix::colSums(cw)
        cws <- pmax(1,Matrix::colSums(cw));
        cw <- t(t(cw)/cws)
        gxd <- Matrix::colSums(cw*arsd$xd)
        gyd <- Matrix::colSums(cw*arsd$yd)
        al <- sqrt(gxd^2+gyd^2);
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
        if(any(vg)) {
          z <- nd %*% cw[,vg]
        } else { NULL }
      },mc.cores=n.cores,mc.preschedule=T))
      cat(". done\n")
      
      
      return(invisible(list(tp=tp,cc=cc,garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,gvel=gv,geshifts=gs,scale=scale)))
    }
    
    
    
    
    } else { # draw individual arrows
      # calculate arrows, draw
      # lapply(1:nrow(emb),function(i) {
      #   # normalized directions to each point
      #   di <- t(t(emb)-emb[i,])
      #   di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
      #   di <- Matrix::colSums(di*tp[,i]) - Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center
      #   
      #   if(fixed.arrow.length) {
      #     suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=0.05,lwd=arrow.lwd))
      #   } else {
      #     ali <- sqrt( (di[1] * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + (di[2]*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2)
      #     suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=min(0.05,ali),lwd=arrow.lwd))
      #   }
      # })
      
      apply(ars,1,function(x) {
        if(fixed.arrow.length) {
          suppressWarnings(arrows(x[1],x[2],x[3],x[4],length=0.05,lwd=arrow.lwd))
        } else {
          ali <- sqrt( ((x[3]-x[1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((x[4]-x[2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2)
          suppressWarnings(arrows(x[1],x[2],x[3],x[4],length=min(0.05,ali),lwd=arrow.lwd))
        }
      })
      
      
    }
  }
  return(invisible(list(tp=tp,cc=cc)))
}



##' Visualize RNA velocities on an existing embedding using Euclidean-based transition probability matrix within the kNN graph.
##'
##'  based on Euclidean distance of the extrapolated cell to others
##' The direction of the arrow is towards n closest neighbors. The magnitude of the arrow is determined by the cosine projection of the velocity on to the chosen direction
##' n=1 will only show arrows for cells that end up being projected closer to some other cell than to the original position
##' n=k (k>1) will show an average direction
##' Given an expression distance between cells d, and ratio of extrapolated to current expression distances between cells f, the transition probability is calculated as exp(- (d*(f^beta))^2/(2*sigma^2) )
##' 
##' @param emb embedding to be used for projection
##' @param vel velocity result
##' @param n neighborhood size (default=30)
##' @param embedding.knn pre-calculated kNN
##' @param cell.colors named color vector for cell plotting
##' @param sigma sigma to use in calculating transition probability from the eucledian distance (estimated automatically by default)
##' @param beta beta parameter used in calculation of transition probability (by default=1)
##' @param arrow.scale additional scaling factor for the arrows (default=1)
##' @param scale scale to use in calculating distances (default: 'log', also supported 'sqrt'
##' @param nPcs number of PCs to project the cells onto (to perform distance calculations in lower dimensions), default=NA which turns off PCA dimensional reduction
##' @param arrow.lwd arrow line width
##' @param xlab x axis label
##' @param ylab y axis label
##' @param control.for.neighborhood.density compensate for cell density variations in the embedding (default: TRUE)
##' @param ntop.trajectories number of top trajectories to trace back for a given cell (when show.trajectories=TRUE)
##' @param do.par whether to reset plotting parameters (default=TRUE)
##' @param show.cell.arrows show detailed velocity projection for the specified cell
##' @param show.cell.trajectories show trajectories for a specified cell
##' @param show.trajectories show top median diffusion trajectories
##' @param show.all.trajectories show all diffusion paths (messy)
##' @param show.cell.diffusion.posterior show diffusion posterior of a given cell
##' @param show.grid.flow show velocity projections on a grid
##' @param diffusion.steps number of diffusion steps to take forward (default=10)
##' @param cell.dist - optional custom distance (must include all of the cells that are intersecting between emb and vel)
##' @param trajectory.spline.shape shape parameter for smoothing median cell trajectories (default=1)
##' @param cell.color.alpha trasparency parameter to apply when showing cell colors
##' @param n.cores number of cores to use in calculations
##' @param n.trajectory.clusters number of trajectory clusters to show median paths for (when show.trajectories=TRUE)
##' @param ... extra parameters are passed to the plot() function
##' @return transition probability matrix
##' @export
show.velocity.on.embedding.eu <- function(emb,vel,n=30,embedding.knn=TRUE,cell.colors=NULL, sigma=NA, beta=1, arrow.scale=1, scale='log', nPcs=NA, arrow.lwd=1, xlab="", ylab="", control.for.neighborhood.density=TRUE, ntop.trajectories=1, do.par=T, show.cell.arrows=NULL, show.cell.trajectories=NULL, show.trajectories=FALSE, show.all.trajectories=FALSE, show.cell.diffusion.posterior=NULL, show.grid.flow=FALSE, diffusion.steps=10, cell.dist=NULL, trajectory.spline.shape=1, cell.color.alpha=0.5, n.cores=defaultNCores(), n.trajectory.clusters=10, ...) {
  em <- vel$current; emn <- vel$projected;
  if(is.null(cell.colors)) { cell.colors <- ac(rep(1,ncol(em)),alpha=0.3); names(cell.colors) <- colnames(em) }
  if(do.par) par(mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  cc <- 'white'
  if(is.null(show.cell.arrows) && is.null(show.cell.diffusion.posterior)) { cc <- cell.colors[rownames(emb)] }
  plot(emb,bg=cc,pch=21,col=ac(1,alpha=cell.color.alpha), xlab=xlab, ylab=ylab, ...);
  
  ccells <- intersect(rownames(emb),colnames(em));
  emn <- emn[,ccells]; em <- em[,ccells]; emb <- emb[ccells,]
  
  if(scale=='log') {
    cat("log scale ... ")
    em <- log10(em+1); emn <- log10(emn+1); 
  } else if(scale=='sqrt') {
    cat("sqrt scale ... ")
    em <- sqrt(em); emn <- sqrt(emn);
  }
  
  if(!is.na(nPcs)) { # run PCA reduction on the em
    cat("reducing to",nPcs,"PCs ... ")
    epc.center <- rowMeans(em);
    epc <- pcaMethods::pca(t(em-epc.center),center=F,nPcs=nPcs);
    em <- t(epc@scores)
    emn <- t(t(emn - epc.center) %*% epc@loadings)
  }
  
  cat("distance ... ")
  cc <- colEuclid(as.matrix(em),as.matrix(emn))
  cc0 <- colEuclid(as.matrix(em),as.matrix(em))
  cd <- (cc0-cc); # reduction in the Euclidean distance

  if(n>nrow(cc)) { n <- nrow(cc) }
  
  # pick reasonable sigma and beta if they weren't provided
  if(is.na(sigma) | is.na(beta)) { mcd <- mean(abs(cd)/cc) }
  # TODO: adaptive methods for signal 
  if(is.na(sigma)) { sigma <- mcd/10 }
  if(is.na(beta)) { beta <- mcd/20 }
  cat("sigma=",round(sigma,3)," beta=",round(beta,3)," transition probs ... ")

  # exp(- (d*(f^beta))^2/(2*sigma^2) )
  f <- (cc/cc0)^beta; diag(f) <- 1;
  tp <- exp(- ((cc0*f)^2) / (2*sigma^2))
  np <- exp(- ((cc0)^2) / (2*sigma^2))

  
  

  if(n<nrow(emb)) {
    if(!is.null(cell.dist)) {
      cat("kNN on provided distance ... ")
      if(!all(labels(cell.dist)==colnames(em))) {
        cat("matching cells between cell.dist and emat/nmat ... ")
        cell.dist <- as.matrix(cell.dist)
        cn <- colnames(em)
        cell.dist <- as.dist(cell.dist[cn,cn]);
      }
      cell.knn <- balancedKNN(t(emb),k=n,maxl=nrow(emb),n.threads=n.cores,dist=cell.dist)
      diag(cell.knn) <- 1;
    } else {
      if(embedding.knn) {
        cat("embedding kNN ... ")
        # define kNNs based on the embedding (L2 distance)
        cell.knn <- balancedKNN(t(emb),k=n,maxl=nrow(emb),dist='euclidean',n.threads=n.cores)
        #diag(cell.knn) <- 0; # disallow self-transitions?
        diag(cell.knn) <- 1;
      } else {
        cat("expression kNN ... ")
        # define kNN based on the correlation distance in high-d
        cell.knn <- balancedKNN(em,k=n,maxl=ncol(em),dist='cor',n.threads=n.cores)
        diag(cell.knn) <- 1;
      }
    }
    tp <- tp*cell.knn;
    np <- np*cell.knn;
  }
  
  # estimate density of the neighborhood
  tp <- t(t(tp)/Matrix::colSums(tp))
  np <- t(t(np)/Matrix::colSums(np))
  #diag(tp) <- diag(np) <- 0;
  #tp.nd <- colSums(tp); np.nd <- colSums(np);
  #tp <- tp * np.nd;

  if(control.for.neighborhood.density) { 
    np.f <- Matrix::diag(np);
    tp <- tp*(np.f)
    np <- np*(np.f)
  }

  #diag(tp) <- diag(np) <- 0;  
  # normalize
  tp <- t(t(tp)/Matrix::colSums(tp))
  np <- t(t(np)/Matrix::colSums(np))

  

  ## if(diffusion.step.size>1) {
  ##   # bring transition probability to the specified power
  ##   require(expm)
  ##   tp <- t(as.matrix(t(tp)) %^% diffusion.step.size);
  ##   np <- t(as.matrix(t(np)) %^% diffusion.step.size);
  ##   tp <- t(t(tp)/Matrix::colSums(tp))
  ##   np <- t(t(np)/Matrix::colSums(np))
  ## }
  
  # normalize transition probabilities
  rownames(tp) <- colnames(tp) <- rownames(np) <- colnames(np) <- colnames(em);
  cat("done\n")
  
  
  if(!is.null(show.cell.diffusion.posterior)) {
    i <- match(show.cell.diffusion.posterior,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # run diffusion
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
    ttp <- t(tp);
    
    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
    }
    cat("done\n");
    # plot
    points(emb,pch=19,col=ac(val2col(cp[show.cell.diffusion.posterior,rownames(emb)],gradient.range.quantile=1),alpha=0.5))
    points(emb[show.cell.diffusion.posterior,1],emb[show.cell.diffusion.posterior,2],pch=3,cex=1,col=1)
  } else if(!is.null(show.cell.arrows)) {
    i <- match(show.cell.arrows,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # plot transition prob for a given cell
    points(emb,pch=19,col=ac(val2col(tp[rownames(emb),show.cell.arrows],gradient.range.quantile=1),alpha=0.5))
    points(emb[i,1],emb[i,2],pch=3,cex=1,col=1)
    di <- t(t(emb)-emb[i,])
    di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
    dir <- Matrix::colSums(di*tp[,i]) 
    dic <- Matrix::colSums(di*np[,i]); # relative to the neighborhood
    dia <- dir-dic;
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],length=0.05,lwd=1,col='blue'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,col='red'))
    suppressWarnings(arrows(emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,lty=1,col='grey50'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dia[1],emb[colnames(em)[i],2]+dia[2],length=0.05,lwd=1,col='black'))
  } else if(show.trajectories) { # show diffusion paths
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
  
    #cpt <- as.array(cp);
    cpl <- list(); cpl[[1]] <- cp;
    #cp <- as.matrix(cp); tp <- as.matrix(tp)
    
    #ep <- as.array(emb)
    ttp <- t(tp);

    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
      #cpt <- abind(cpt,cp,along=3)
      cpl[[i+1]] <- cp;
      # clean up to zero out all but top n cells
      #cp <- t(apply(cp,1,function(x) { x[x<sort(x,decreasing=TRUE)[10]] <- 0; x }))
      #diag(cp) <- 0; #  prohibit the cell from returning to itself
      #cp <- cp/Matrix::rowSums(cp)
    }


    #cpt <- abind(lapply(cpl,as.matrix),along=3)
    
    # calculate probabilistic trajectories to the final ntop points
    
    # rank final points by probability
    cpo <- t(apply(-cp,1,order))

    # graph-based walkback approach
    
    # construct a walkback graph
    trp <- as(ttp,'dgTMatrix')
    cat("constructing path graph ... ")
    x <- do.call(rbind,lapply(1:(diffusion.steps+1),function(i) {
      cbind(i=trp@i+1 + (i-1)*nrow(cp), # current time step
            j=trp@j+1 + (i)*nrow(cp))
    }))
    x <- x[x[,2]<=nrow(cp)*(diffusion.steps+1),]
    x <- spMatrix(nrow(cp)*(diffusion.steps+1),nrow(cp)*(diffusion.steps+1),i=x[,1],j=x[,2],x=rep(1,nrow(x)))
    g <- igraph::graph.adjacency(x,mode='directed')
    rm(x); gc();
    
    # find topn trajectories for each cell
    cat("tracing shortest trajectories ... ")
    sps <- parallel::mclapply(1:nrow(cp),function(celli) {
      top.desti <- order(cp[celli,],decreasing=TRUE)[1:ntop.trajectories]
      # calculate cell-specific weights
      cw <- unlist(lapply(cpl,function(d) as.numeric(trp@x*(d[celli,trp@i+1]))))
      cw <- cw[1:igraph::ecount(g)] # trim extra edges
      # convert into penalty scores
      cw <- -log(cw)
      sp <- igraph::shortest_paths(g,from=celli,to=nrow(cp)*(diffusion.steps-1)+top.desti,weights=cw,mode='out')
      # remove time offset on the path nodes
      sp <- lapply(sp$vpath,function(x) { y <- (as.integer(x) %% nrow(cp)); y[y==0] <- nrow(cp); y});
      names(sp) <- rownames(cp)[top.desti]
      sp
    },mc.cores=n.cores)



    # cluster paths
    cat("clustering ... ")
    all.cells <- 1:nrow(cp)
    #spuci <- do.call(cbind,lapply(sps,function(x) all.cells %in% x[[1]]))
    # filter out empty paths
    sps <- lapply(sps,function(y) y[unlist(lapply(y,function(x) length(unique(x))))>1])
    spuci <- do.call(cbind,lapply(sps,function(y) do.call(cbind,lapply(y,function(x) all.cells %in% x))))
    usps <- unlist(sps,recursive=F); # will be used in looking up median trajectories in plotting

    
    spuci.dist <- as.matrix(dist(t(spuci),method = 'manhattan'))
    spuci.pam <- pam(spuci.dist,n.trajectory.clusters)
    cat("done.\n")
    
    # bezier
    # determine common start/end points
    #plot(emb,bg='white',pch=21,col=ac(1,alpha=cell.color.alpha), xlab=xlab, ylab=ylab);
    lapply(1:length(spuci.pam$id.med),function(cn) {
      if(length(usps[[spuci.pam$id.med[cn]]])>0){
        mp <- usps[[spuci.pam$id.med[cn]]]; mp <- mp[!duplicated(mp)]
        bp <- data.frame(do.call(cbind,xspline(emb[mp,],shape=trajectory.spline.shape,draw=F)))
        lines(bp$x,bp$y,col=ac(1,alpha=0.6))
        bp <- bp[abs(diff(bp$x))+abs(diff(bp$y))>1e-5,]
        ai <- round(length(bp$x)*c(0.2,0.8,0.5))
        arrows(bp$x[ai],bp$y[ai],bp$x[ai+1],bp$y[ai+1],angle=30,length=0.1,col=ac(1,alpha=0.6))
      }
    })
    return(invisible(list(spuci.dist=spuci.dist,sps=sps,tp=tp,cpl=cpl)))
  } else if(!is.null(show.cell.trajectories)) {
    # show optimal path(s) for a particular cell
    celli <- match(show.cell.trajectories,rownames(emb));
    if(is.na(celli)) stop(paste('specified cell',show.cell.trajectories,'is not in the embedding'))

    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities
    rownames(cp) <- colnames(cp) <- rownames(tp);
    
    #cpt <- as.array(cp);
    cpl <- list(); cpl[[1]] <- cp;
    #cp <- as.matrix(cp); tp <- as.matrix(tp)
    
    #ep <- as.array(emb)
    ttp <- t(tp);
    
    # run diffusion steps to figure out end positions
    cat("simulating diffusion ... ")
    for(i in 1:diffusion.steps) {
      cp <- cp %*% ttp;
      #cp[cp<1e-5] <- 0;
      cpl[[i+1]] <- cp;
    }

    # rank final points by probability
    cpo <- t(apply(-cp,1,order))

    # graph-based walkback approach
    
    # construct a walkback graph
    trp <- as(ttp,'dgTMatrix')
    
    cat("constructing path graph ... ")
    x <- do.call(rbind,lapply(1:(diffusion.steps+1),function(i) {
      cbind(i=trp@i+1 + (i-1)*nrow(cp), # current time step
            j=trp@j+1 + (i)*nrow(cp))
    }))
    x <- x[x[,2]<=nrow(cp)*(diffusion.steps+1),]
    x <- spMatrix(nrow(cp)*(diffusion.steps+1),nrow(cp)*(diffusion.steps+1),i=x[,1],j=x[,2],x=rep(1,nrow(x)))
    g <- igraph::graph.adjacency(x,mode='directed')
    rm(x); gc();
    
    # find topn trajectories for each cell
    cat("tracing shortest trajectories ... ")
    top.desti <- order(cp[celli,],decreasing=TRUE)[1:ntop.trajectories]
    
    
    # calculate cell-specific weights
    cw <- unlist(lapply(cpl,function(d) as.numeric(trp@x*(d[celli,trp@i+1]))))
    cw <- cw[1:igraph::ecount(g)] # trim extra edges
    # convert into penalty scores
    cw <- -log(cw)
    sp <- igraph::shortest_paths(g,from=celli,to=nrow(cp)*(diffusion.steps)+top.desti,weights=cw,mode='out')
    # remove time offset on the path nodes
    sp <- lapply(sp$vpath,function(x) { y <- (as.integer(x) %% nrow(cp)); y[y==0] <- nrow(cp); y});
    names(sp) <- rownames(cp)[top.desti]
    cat("done.\n")
    lapply(sp,function(mp) {
      if(!is.null(mp) && length(mp)>0) {
        mp <- mp[!duplicated(mp)]
        if(length(mp)>1)  {
          lines(emb[mp,1],emb[mp,2],col=8,lty=3)
          bp <- data.frame(do.call(cbind,xspline(emb[mp,],shape=trajectory.spline.shape,draw=F)))
          lines(bp$x,bp$y,col=ac(1,alpha=0.6))
          bp <- bp[abs(diff(bp$x))+abs(diff(bp$y))>1e-5,]
          ai <- round(length(bp$x)*c(0.2,0.8,0.5))
          arrows(bp$x[ai],bp$y[ai],bp$x[ai+1],bp$y[ai+1],angle=30,length=0.1,col=ac(1,alpha=0.6))
        }
      }
    })
    return(invisible(sp))
  } else if(show.all.trajectories) { # show diffusion paths
    cp <- Matrix::Diagonal(ncol(tp)); # cell position probabilities row-from col-to
    rownames(cp) <- colnames(cp) <- rownames(tp);
    ep <- as.array(emb)
    
    for(i in 1:diffusion.steps) {
      cp <- cp %*% t(tp);
      # expected position
      cpm <- t(apply(cp,1,function(x) { x[x<sort(x,decreasing=TRUE)[3]] <- 0; x/sum(x) }))
      epi <- as.matrix(cpm %*% emb);
      #epi <- as.matrix(cp %*% emb);
      ep <- abind::abind(ep,epi,along=3)
    }
    apply(ep,c(1),function(d) {
      lines(d[1,],d[2,],lwd=1,col=ac(1,alpha=0.05))
    })
  } else {
    # calculate arrows, draw
    lapply(1:nrow(emb),function(i) {
      # normalized directions to each point
      di <- t(t(emb)-emb[i,])
      di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
      di <- Matrix::colSums(di*tp[,i]) - Matrix::colSums(di*np[,i]); # relative to expected kNN center
      suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=0.05,lwd=arrow.lwd))
    })
  }
  return(invisible(tp))
  
}


##' DEPRECATED: Read in cell-specific bam files for SMART-seq2 measurement
##' This function is deprecated. Please use velocyto.py to prepare loom file from SMART-seq2 bam files.
##' @title read.smartseq2.bams
##' @param bam.files list of bam files
##' @param annotation.file refFlat genome annotation file (use gtfToGenePred to generate refFlat file from gtf)
##' @param min.exon.count minimum number of reads (across all cells) for an exon to be considered expressed in the dataset
##' @param n.cores number of cores to use
##' @return a list containing: emat - exonic (spliced) read count matrix ; iomat - intronic (unspliced) matrix; smat - spanning read matrix; base.df - data frame containing gene structural information; exons - exon annotation and read counts; genes - gene annotation table with additional structural info; expr.lstat - gene length statistics when considering only expressed exons
##' @export
read.smartseq2.bams <- function(bam.files,annotation.file,min.exon.count=100,n.cores=defaultNCores()) {
  # read in annotation 
  # TODO: enable direct gtf read
  cat("reading gene annotation ... ")
  x <- read.delim(annotation.file,header=F,sep="\t",stringsAsFactors=F)
  genes <- data.frame(name=x[,1],chr=x[,3],strand=x[,4],start=x[,5],end=x[,6],stringsAsFactors=F)
  genes$p5 <- genes$start; genes$p5[genes$strand=="-"] <- genes$end[genes$strand=="-"];
  genes$p3 <- genes$end; genes$p3[genes$strand=="-"] <- genes$start[genes$strand=="-"];
  genes$size <- genes$end - genes$start;
  cat("done (",nrow(genes),"genes)\n")
  
  
  # parse out exon information
  cat("parsing exon information ... ")
  exons <- do.call(rbind,lapply(1:nrow(x),function(i) {
    df <- do.call(cbind,strsplit(as.character(x[i,c(10,11)]),','))
    cbind(df,rep(as.character(x[i,1]),nrow(df)))
  }))
  exons <- data.frame(gene=as.character(exons[,3]),start=as.numeric(exons[,1]),end=as.numeric(exons[,2]),stringsAsFactors=F)
  exons$chr <- genes$chr[match(exons$gene,genes$name)]
  
  # eliminate duplicated exons? - points.within will count reads only once anyhow
  exons <- exons[!duplicated(paste(exons[,1],exons[,2],exons[,3])),]
  
  # eliminate multiple variants - keep the largest ones
  genes <- genes[order(genes$size,decreasing=T),]
  genes <- genes[!duplicated(genes$name),]
  cat("done\n")
  
  # read in expression data
  #  - ultimately we'll get three kinds of matrices here:
  #      1. exonic reads
  #      2. intronic-only reads
  #      3. spanning reads (overlapping an intron and an exon)
  #  - we will also use bam files to count the number of expressed exons,
  #    to adjust the structural parameters of each gene. Though this is not
  #    as important for the intron-only models.
  
  # read in all bam files
  cat("reading in",length(bam.files),"bam files ... ")
  # annotate individual reads
  cdl <- parallel::mclapply(bam.files,t.annotate.bam.reads,genes=genes,exons=exons,margin=1,exon.margin=1,mc.cores=n.cores)
  cat("done\n")
  # get count estimates per gene
  cat("estimating gene counts ... ")
  edl <- parallel::mclapply(cdl,t.get.estimates2,genes=genes,mc.cores=n.cores)
  cat("done\n")

  cat("adjusting gene annotation based on expressed regions ... ")
  # count total number of reads per exon to get the ones that are expressed ...
  tl <- Matrix::colSums(do.call(rbind,parallel::mclapply(cdl,function(x) { ect <- table(c(x$exonstart,x$exonend)); fect <- rep(0,nrow(exons)); fect[as.integer(names(ect))] <- ect; fect; },mc.cores=n.cores,mc.preschedule=T)))
  
  # calculate gene length statistics, based on the expressed exons
  expr.exons <- exons[tl>min.exon.count,]; # expressed exons - those showing >100 reads dataset-wide
  expr.lstat <- lengthstats2(1e8,genes=genes,exons=expr.exons); # intron/exon length statistics for genes, considering only expressed exons
  
  # compile joint table
  df <- data.frame(il=log10(expr.lstat[,2]+1),el=log10(expr.lstat[,3]+1)); rownames(df) <- rownames(expr.lstat)
  df$nex <- as.integer(table(expr.exons$gene)[rownames(df)]);
  df$nex[is.na(df$nex)] <- 0
  cat("done\n")
  
  # construct sparse input matrices
  # exonic counts
  emat <- do.call(cbind,lapply(edl,function(d) { (d[,'exon']) }));  emat[!is.finite(emat)] <- 0;  emat <- as(emat,'dgCMatrix')
  # spanning counts
  smat <- do.call(cbind,lapply(edl,function(d) { (d[,'span']) }));  smat[!is.finite(smat)] <- 0;  smat <- as(smat,'dgCMatrix')
  # intron-only counts
  iomat <- do.call(cbind,lapply(edl,function(d) { (d[,'introno']) }));  iomat[!is.finite(iomat)] <- 0;  iomat <- as(iomat,'dgCMatrix')
  
  
  return(list(emat=emat,iomat=iomat,smat=smat,base.df=df,exons=exons,genes=genes,expr.lstat=expr.lstat))
}

##' Read in cell-specific bam files for STRT/C1
##'
##' @param bam.files list of bam files (one per cell)
##' @param annotation.file gene annotation refFlat file
##' @param min.exon.count minimum number of molecules (across all cells) for the exon to be considered expressed
##' @param n.cores number of cores to use
##' @param min.umi.reads minimum number of read required per UMI/gene combination to be counted (defaults to 1)
##' @return a list structure analogous to the return of read.smartseq2.bams(), counting molecules instead of reads.
read.strtc1.bams <- function(bam.files,annotation.file,min.exon.count=100,n.cores=defaultNCores(),min.umi.reads=1) {
  # read in annotation 
  # TODO: enable direct gtf read
  cat("reading gene annotation ... ")
  x <- read.delim(annotation.file,header=F,sep="\t",stringsAsFactors=F)
  genes <- data.frame(name=x[,1],chr=x[,3],strand=x[,4],start=x[,5],end=x[,6],stringsAsFactors=F)
  genes$p5 <- genes$start; genes$p5[genes$strand=="-"] <- genes$end[genes$strand=="-"];
  genes$p3 <- genes$end; genes$p3[genes$strand=="-"] <- genes$start[genes$strand=="-"];
  genes$size <- genes$end - genes$start;
  cat("done (",nrow(genes),"genes)\n")
  
  
  # parse out exon information
  cat("parsing exon information ... ")
  exons <- do.call(rbind,lapply(1:nrow(x),function(i) {
    df <- do.call(cbind,strsplit(as.character(x[i,c(10,11)]),','))
    cbind(df,rep(as.character(x[i,1]),nrow(df)))
  }))
  exons <- data.frame(gene=as.character(exons[,3]),start=as.numeric(exons[,1]),end=as.numeric(exons[,2]),stringsAsFactors=F)
  exons$chr <- genes$chr[match(exons$gene,genes$name)]
  
  # eliminate duplicated exons? - points.within will count reads only once anyhow
  exons <- exons[!duplicated(paste(exons[,1],exons[,2],exons[,3])),]
  
  # eliminate multiple variants - keep the largest ones
  genes <- genes[order(genes$size,decreasing=T),]
  genes <- genes[!duplicated(genes$name),]
  cat("done\n")
  
  # read in expression data
  #  - ultimately we'll get three kinds of matrices here:
  #      1. exonic reads
  #      2. intronic-only reads
  #      3. spanning reads (overlapping an intron and an exon)
  #  - we will also use bam files to count the number of expressed exons,
  #    to adjust the structural parameters of each gene. Though this is not
  #    as important for the intron-only models.
  
  # read in all bam files
  cat("reading in",length(bam.files),"bam files ... ")
  # annotate individual reads

  #t.reduce.umi <- function(z) { z[!duplicated(paste(gsub(".*?_","",z$name),z$gene)),] }
  t.reduce.umi <- function(z) {
    z$umig <- paste(gsub(".*?_","",z$name),z$gene)
    z[match(names(which(table(z$umig)>=min.umi.reads)),z$umig),]
  }
  
  cdl <- parallel::mclapply(bam.files,function(x) t.reduce.umi(t.annotate.bam.reads(x,genes=genes,exons=exons,margin=1,exon.margin=1,use.names=T)),mc.cores=n.cores)
  
  
  cat("done\n")
  # get count estimates per gene
  cat("estimating gene counts ... ")
  edl <- parallel::mclapply(cdl,t.get.estimates2,genes=genes,mc.cores=n.cores)
  cat("done\n")

  cat("adjusting gene annotation based on expressed regions ... ")
  # count total number of reads per exon to get the ones that are expressed ...
  tl <- Matrix::colSums(do.call(rbind,parallel::mclapply(cdl,function(x) { ect <- table(c(x$exonstart,x$exonend)); fect <- rep(0,nrow(exons)); fect[as.integer(names(ect))] <- ect; fect; },mc.cores=n.cores,mc.preschedule=T)))
  
  # calculate gene length statistics, based on the expressed exons
  expr.exons <- exons[tl>min.exon.count,]; # expressed exons - those showing >100 reads dataset-wide
  expr.lstat <- lengthstats2(1e8,genes=genes,exons=expr.exons); # intron/exon length statistics for genes, considering only expressed exons
  
  # compile joint table
  df <- data.frame(il=log10(expr.lstat[,2]+1),el=log10(expr.lstat[,3]+1)); rownames(df) <- rownames(expr.lstat)
  df$nex <- as.integer(table(expr.exons$gene)[rownames(df)]);
  df$nex[is.na(df$nex)] <- 0
  cat("done\n")
  
  # construct sparse input matrices
  # exonic counts
  emat <- do.call(cbind,lapply(edl,function(d) { (d[,'exon']) }));  emat[!is.finite(emat)] <- 0;  emat <- as(emat,'dgCMatrix')
  # spanning counts
  smat <- do.call(cbind,lapply(edl,function(d) { (d[,'span']) }));  smat[!is.finite(smat)] <- 0;  smat <- as(smat,'dgCMatrix')
  # intron-only counts
  iomat <- do.call(cbind,lapply(edl,function(d) { (d[,'introno']) }));  iomat[!is.finite(iomat)] <- 0;  iomat <- as(iomat,'dgCMatrix')
  
  
  return(list(emat=emat,iomat=iomat,smat=smat,base.df=df,exons=exons,genes=genes,expr.lstat=expr.lstat))
}



# parse bam file and annotate the reads relative to genes/exons
t.annotate.bam.reads <- function(fname, genes, exons, chrl=unique(genes$chr), test.strand=F, margin=3e3, tags=NULL,use.names=FALSE,exon.margin=0) {

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package \"GenomicRanges\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package \"Rsamtools\" needed for this function to work. Please install it.",call. = FALSE)
  }

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package \"GenomeInfoDb\" needed for this function to work. Please install it.",call. = FALSE)
  }

  
  if(!is.null(tags)) {
    param <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),tag=unlist(tags))
  } else {
    param <- Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE))
  }
  z <- GenomicAlignments::readGAlignments(fname,param=param,use.names=use.names)
  
  bam.data <- data.frame(chr=as.vector(GenomeInfoDb::seqnames(z)),start=BiocGenerics::start(z),end=BiocGenerics::end(z),strand=as.vector(BiocGenerics::strand(z)),stringsAsFactors=F)
  ## if(!is.null(tags)) {
  ##   bam.data <- cbind(bam.data,as.data.frame(S4Vectors::elementMetadata((z))))
  ## }
  if(use.names) {
    bam.data$name <- names(z)
  }
  bam.data <- bam.data[bam.data$chr %in% chrl,]
  
  chrl <- chrl[chrl %in% unique(bam.data$chr)]
  
  ## assign reads to genes
  x <- do.call(rbind,lapply(chrl,function(chr) {
    ri <- which(bam.data$chr==chr)
    results <- data.frame(bam.data[ri,],gene=rep(NA,length(ri)), exonstart=0, exonend=0, readtype='NI',stringsAsFactors=F) ## intronic read unless evidence to change it
    gi <- which(genes$chr==chr)
    ei <- which(exons$chr==chr)
    
    # check if the read is associated with any gene
    pi1 <- points.within(bam.data$start[ri],genes$start[gi]-margin,genes$end[gi]+margin); # note that many gene fragments are supplied at once - no need to loop in R!
    pi2 <- points.within(bam.data$end[ri],genes$start[gi]-margin,genes$end[gi]+margin);
    vi <- pi1>0 | pi2>0;
    if(any(vi)) { 
      results$gene[vi] <- genes$name[gi[pmax(pi1[vi],pi2[vi])]]; # take largest matched gene index .. should check for consistency
    }
    
    # check exon mapping
    pi1 <- points.within(bam.data$start[ri],exons$start[ei]-exon.margin,exons$end[ei]+exon.margin);
    pi2 <- points.within(bam.data$end[ri],exons$start[ei]-exon.margin,exons$end[ei]+exon.margin);
    
    # see if both ends map to the same exon - that's NC
    vi <- pi1>0  & pi1==pi2
    if(any(vi)) { results$readtype[vi] <- 'NC'; results$gene[vi] <- exons$gene[ei[pi1[vi]]]; }
    
    # different exons - NE
    vi <- pi1>0 & pi2>0 & pi1!=pi2
    if(any(vi)) { results$readtype[vi] <- 'NE'; results$gene[vi] <- exons$gene[ei[pi1[vi]]]; }
    
    # one exon, and one something else - assume NS for now, will check for margins later
    vi <- sign(pi1) != sign(pi2);
    if(any(vi)) { results$readtype[vi] <- 'NS'; results$gene[vi] <- exons$gene[ei[pmax(pi1,pi2)[vi]]]; }
    
    # note: I am not sure what these two values are needed for, but here's how to get them
    pi1[pi1==-1] <- NA; pi2[pi2==-1] <- NA;
    results$exonstart <- ei[pi1]; results$exonend <- ei[pi2];
    
    # assign margin classes (N5 - fully in the 5' margin, N3 - fully in the 3', N5b - spanning exon and 5' margin, N3b - spanning exon and 3' margin)
    if(margin>0) {
      # start is in the start margin
      pi1 <- points.within(bam.data$start[ri],genes$start[gi]-margin,genes$start[gi]-1);
      pi1check <- points.within(bam.data$start[ri],genes$start[gi],genes$end[gi]); # check for genes
      ivi <- pi1>0 & pi1check>0 # both in margin and in some other gene - should be invalidated
      if(any(ivi)) { results$gene[ivi] <- NA }
      
      # if it is in vi, then it must be partially margin (joining margin and the first/last exon), otherwise it's fully in margin
      vi2 <- pi1>0 & vi;
      if(any(vi2)) { results$readtype[vi2] <- ifelse(genes$strand[gi[pi1[vi2]]]=='+','N5b','N3b') }
      vi2 <- pi1>0 & !vi;
      if(any(vi2)) { results$readtype[vi2] <- ifelse(genes$strand[gi[pi1[vi2]]]=='+','N5','N3') }
      
      # and now check the other end of the read
      # end is in the margin
      pi2 <- points.within(bam.data$end[ri],genes$end[gi]+1,genes$end[gi]+margin);
      pi2check <- points.within(bam.data$end[ri],genes$start[gi],genes$end[gi]); # check for genes
      ivi <- pi2>0 & pi2check>0 # both in margin and in some other gene - should be invalidated
      if(any(ivi)) { results$gene[ivi] <- NA }
      
      # if it is in vi, then it must be partially margin (joining margin and the first/last exon), otherwise it's fully in margin
      vi2 <- pi2>0 & vi;
      if(any(vi2)) { results$readtype[vi2] <- ifelse(genes$strand[gi[pi2[vi2]]]=='+','N3b','N5b') }
      vi2 <- pi2>0 & !vi;
      if(any(vi2)) { results$readtype[vi2] <- ifelse(genes$strand[gi[pi2[vi2]]]=='+','N3','N5') }
    }
    
    # omit reads that didn't get assinged to genes (or were invalidated)
    results <- results[!is.na(results$gene),]
    return(results)
  }))
}

# count different read types per gene
t.get.estimates2 <- function(cd,genes) {
  cd$gene <- match(cd$gene,genes$name);
  # number of reads within the exons
  vi <- cd$readtype %in% c("NC","NE");  exon.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  
  # intronic reads
  vi <- cd$readtype %in% c("NI","NS"); intron.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  # intron only counts
  vi <- cd$readtype %in% c("NI"); introno.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  # span only counts
  vi <- cd$readtype %in% c("NS"); span.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  
  # number of reads in the upstream region (don't include those briding exon, since they're supposed to be a different strand)
  vi <- cd$readtype %in% c("N5"); upstream.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  # downstream region (running from the exon is fine, since it can be the same strand)
  vi <- cd$readtype %in% c("N3","N3b"); downstream.counts <- tapply(cd$gene[vi],cd$gene[vi],length)
  
  
  # construct report matrix
  df <- matrix(NA,nrow=nrow(genes),ncol=7)
  df[as.integer(names(exon.counts)),1] <- exon.counts
  df[,2] <- 0
  # normalize by the expected fraction of reads from the last 1kb (assume that the last exon is ~1kb) TODO: transcript-coordinate estimation
  #df[,2] <- df[,2]/lstat.e3$e*lstat$e;
  
  df[as.integer(names(intron.counts)),3] <- intron.counts
  df[as.integer(names(upstream.counts)),4] <- upstream.counts
  df[as.integer(names(downstream.counts)),5] <- downstream.counts
  df[as.integer(names(introno.counts)),6] <- introno.counts
  df[as.integer(names(span.counts)),7] <- span.counts
  rownames(df) <- genes$name; colnames(df) <- c("exon","mg5","intron","upstream","downstream","introno","span")
  df
}


lengthstats2 <- function(maxlen,genes,exons,n.cores=defaultNCores(),p3=FALSE) {
  lf <- do.call(rbind,parallel::mclapply(1:nrow(genes),function(gi) {
    ei <- which(exons$gene==genes$name[gi])
    if((genes$strand[gi]=="+" & !p3) | (genes$strand[gi]=="-" & p3)) {
      nt <- seq(genes$start[gi],min(genes$start[gi]+maxlen-1,genes$end[gi]),by=1);
    } else {
      nt <- seq(max(genes$end[gi]-maxlen-1,genes$start[gi]),genes$end[gi],by=1)
    }
    c(t=length(nt),i=sum(points.within(nt,exons$start[ei],exons$end[ei],sorted=T)==-1))
  },mc.cores=n.cores,mc.preschedule=T))
  lf <- cbind(lf,e=lf[,"t"]-lf[,"i"])
  rownames(lf) <- genes$name
  data.frame(lf)
}

# quick utility function to calculate library sizes using edgeR
edgeR.libsize <- function(mat, ...) {
  f <- edgeR::calcNormFactors(mat)
  f <- f/exp(mean(log(f)))
  Matrix::colSums(mat)*f
}

# quick self-naming vector routine
sn <- function(x) { names(x) <- x; x}

#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...); names(y) <- names(x); return(y)}

# quick function to map value vector to colors
val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(x)>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(na.omit(x),p=c(1-gradient.range.quantile,gradient.range.quantile)))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(na.omit(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(na.omit(abs(x)),p=gradient.range.quantile))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(na.omit(max(abs(x))))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  }
  
  gp <- gradientPalette[x*(length(gradientPalette)-1)+1]
  if(!is.null(names(x))) { names(gp) <- names(x) }
  gp
}

# estimate expression of a projected cell given the original expression matrix and deltaE matrix
# delta is the amount of time onto which the projection should be done
t.get.projected.cell <- function(em,deltae,delta=1,target.mult=1e3,model.mult=1e3,size.normalize=FALSE) {
  rz <- matrix(0,nrow=nrow(em),ncol=ncol(em)); colnames(rz) <- colnames(em); rownames(rz) <- rownames(em)
  gn <- intersect(rownames(deltae),rownames(rz))
  rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,]*target.mult/model.mult; # correcting for different mult
  emn <- em+delta*rz;  emn[emn<0] <- 0
  # rescale ... note: none of these seem appropriate
  #i <- 11;smoothScatter(sqrt(em[,i]),sqrt(emn[,i])); abline(a=0,b=1)
  #emn.scale <- unlist(lapply(1:ncol(em),function(i) { n <- calcNormFactors(cbind(em[,1],emn[,1])); n[2]/n[1]}))
  #emn <- t(t(emn)/emn.scale)
  if(size.normalize) {
    emn <- t(t(emn)/Matrix::colSums(emn)*Matrix::colSums(em))
  }
  emn
}

# calculates the difference in the number of counts based on the library size, renormalizes
# note: also introduces centripetal velocity
t.get.projected.cell2 <- function(em,cellSize,deltae,mult=1e3,delta=1) {
  rz <- matrix(0,nrow=nrow(em),ncol=ncol(em)); colnames(rz) <- colnames(em); rownames(rz) <- rownames(em)
  gn <- intersect(rownames(deltae),rownames(rz))
  rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,]; 
  # translate fpm delta into the number of molecules based on the current cell size
  rz <- t(t(rz)*cellSize)
  emm <- t(t(em)*cellSize)
  emn <- emm + rz*delta; 
  emn[emn<0] <- 0;
  newCellSize <- (cellSize+Matrix::colSums(emn-emm)/mult)
  emn <- t(t(emn)/newCellSize)
  
  #emn <- t(t(emn)/Matrix::colSums(emn)*Matrix::colSums(em))
  emn
}

# reads layers/spliced/unspliced/ambiguous from a loom file
##' Read in loom matrices of spliced/unpsliced reads as
##' prepared by velocyto.py CLI
##'
##' @param file loom file name
##' @param engine Use hdf5r or h5 to import loom file
##' @return a list containing spliced, unspliced, ambiguous and spanning matrices
##' @export
read.loom.matrices <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}


##' identify positions of likely internal priming sites by looking for polyA/polyT stretches within
##' annotated intronic regions
##'
##' @param gtf.file location of a gtf file with gene annotations
##' @param genome bioC genome structure (i.e. Mmusculus)
##' @param genome.name name of the genome assembly (i.e. mm10)
##' @param w A/T weight in the PWM
##' @param n length of the motif
##' @param min.score minimal required match score
##' @param add.chr whether to add 'chr' prefix to the chromosome names in the gtf annotation (to match bioC)
##' @return a data frame containing list of likely internal priming sites, listing chromosome ($chr), positions ($start/$end), name, PWM matching score, PWM match strand ($strand), gene, gene strand ($gs), and whether the motif is in concordant direction of gene transcription ($conc)
##' @examples
##' \dontrun{
##' library(BSgenome.Mmusculus.UCSC.mm10)
##' ip.mm10 <- find.ip.sites('refdata-cellranger-mm10-1.2.0/genes/genes.gtf',Mmusculus,'mm10')
##' }
##' @export
find.ip.sites <- function(gtf.file,genome,genome.name,w=0.9,n=15,min.score='80%',add.chr=TRUE) {
  
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package \"GenomicRanges\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package \"Biostrings\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # helper functions to move out of bioC classes
  grange2df <- function(x,name='hit') {
    data.frame(chr=as.character(GenomicRanges::seqnames(x)),
               start=GenomicRanges::start(x),
               end=GenomicRanges::end(x),
               name=name,
               score=GenomicRanges::score(x),
               strand=GenomicRanges::strand(x))
  }
  
  # read specified features from gtf file, recording specified attributes
  t.read.gtf <- function(file,feature="gene",atts=c("gene_id"),n.cores=30) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package \"data.table\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    #x <- read.table(file,header=F,sep="\t",stringsAsFactors=F)
    x <- data.table::fread(file,header=F,sep="\t",stringsAsFactors=F,data.table=FALSE)
    vi <- which(x[,3]==feature)
    names(atts) <- atts;
    ad <- do.call(cbind,lapply(atts,function(att) {
      gsub("\\\"","",gsub(paste('.*',att,' ([^;]*);.*',sep=""),"\\1",x[vi,9]))
    }))
    df <- x[vi,c(1,4,5,7)]
    colnames(df) <- c("chr","start","end","strand");
    cbind(df,ad,stringsAsFactors=F)
  }

  cat("reading genes ... ");
  genes <- t.read.gtf(gtf.file,atts=c("gene_name"))
  cat("done\n");
  cat("reading exons ... ");
  exons <- t.read.gtf(gtf.file,feature='exon',atts=c("gene_name","transcript_biotype"))
  exons <- exons[exons$transcript_biotype=='protein_coding',] # want to mask everything else
  if(add.chr) {
    genes$chr <- paste0('chr',genes$chr);
    exons$chr <- paste0('chr',exons$chr);
  }
  cat("done\n");
  cat("making ranges ... ")
  gene.ranges <- GenomicRanges::GRanges(genes$chr,IRanges::IRanges(start=genes[,2],end=genes[,3]),strand=genes$strand)
  names(gene.ranges) <- genes$gene_name
  exon.ranges <- GenomicRanges::GRanges(exons$chr,IRanges::IRanges(start=exons[,2],end=exons[,3]),strand=exons$strand)
  cat("done\n");
  cat("matching hits ... ")
  N <- 1e3;
  m <- as.matrix(rbind(A=rep(round(w*N),n),C=rep(round((1-w)/3*N),n),G=rep(round((1-w)/3*N),n),T=rep(round((1-w)/3*N),n))); storage.mode(m) <- 'integer'
  pwm <- Biostrings::PWM(m)
  hits <- Biostrings::matchPWM(pwm,genome,min.score=min.score,with.score=T)
  cat("done\n");
  cat("annotating hits .");
  df <- grange2df(hits,name=paste0('(A)',n)) # convert into a table
  cat(".");
  df <- do.call(rbind,tapply(1:nrow(df),df$chr,function(ii) df[ii[order(df$start[ii])],])) # sort
  cat(".");
  df <- rbind(groupMotifs(df[df$strand=='+',]),groupMotifs(df[df$strand=='-',])) # cluster
  cat(".");
  df <- do.call(rbind,tapply(1:nrow(df),df$chr,function(ii) df[ii[order(df$start[ii])],])) # sort again
  cat(".");
  # get gene strand information and classify concordance
  sn <- function(x) { names(x) <- x; x}
  gene.margin <- 0;
  df <- do.call(rbind,lapply(sn(intersect(unique(genes$chr),unique(df$chr))),function(chr) {
    vg <- which(genes$chr==chr);
    vm <- which(df$chr==chr);
    ei <- points.within((df$start[vm]+df$end[vm])/2,genes$start[vg]-gene.margin,genes$end[vg]+gene.margin)
    vi <- which(ei>-1)
    x <- cbind(df[vm[vi],],gene=genes$gene_name[vg[ei[vi]]],gs=genes$strand[vg[ei[vi]]])
    x$conc <- x$strand==x$gs
    x
  }))
  cat(". done\n");
  return(df)
}

##' read in detailed molecular mapping info from hdf5 file as written out by "-d" option of velocyto.py
##'
##' @param fname name of the hdf5 detailed molecular mapping (debug) file, written out by velocyto.py
##' @param cell.clusters optional cell cluster factor
##' @param internal.priming.info optionall internal priming info, as produced by find.ip.sites() function
##' @param min.exon.count minimal total (dataset-wide) number of molecules for an exon to be considered expressed
##' @param n.cores number of cores to use
##' @param engine use either h5 or hdf5r to read the input hdf5 file
##' @return a list containing gene structural information data structure ($gene.df, with el,il, nex,nipconc,nipdisc columns corresponding to the log10 exonic length, intronic length, number of exons, numebr of internal concordant and discordant priming sites, respectively), and $info tables from the hdf5 file with an additional per-cluster entry $cluster.feature.counts table showing per-feature (rows) per-cluster (column) molecule counts (if cell.clusters are not supplied $info$cluster.feauture.counts will contain one column, 'all' giving dataset-wide counts)
##' @export
read.gene.mapping.info <- function(fname,cell.clusters=NULL,internal.priming.info=NULL,min.exon.count=10,n.cores=defaultNCores(),engine='hdf5r') {
  cat("reading in mapping info from",fname,' via', engine)
  if (engine == 'h5') {
    f <- h5::h5file(fname,mode='r')
  } else if (engine == 'hdf5r') {
    f <- hdf5r::H5File$new(fname, mode='r')
  } else {
    stop('Unknown engine. Use either hdf5r or h5.\n')
  }
  #list.datasets(f)
  # read in info tables
  if (engine == 'h5') info <- lapply(sn(c("chrm","exino","features_gene","is_intron","is_last3prime","start_end","strandplus","tr_id")),function(n) { cat('.'); f[paste('/info',n,sep='/')][] })
  if (engine == 'hdf5r') {
  info <- lapply(
    sn(c("chrm","exino","features_gene", "is_intron","is_last3prime",
         "start_end","strandplus","tr_id")),
    function(n) {
      cat('.'); x <- f[[paste('info', n, sep='/')]]$read();
      if (is.matrix(x)) {t(x)} else {x} })
  }
  info$chrm <- gsub("^chr","",info$chrm)
  cat(" done\n")
  # extract cell names
  if (engine == 'h5') cnames <- gsub('/pos','',gsub('/cells/','',grep('/pos',grep("/cells/", h5::list.datasets(f),value=T),value=T)))
  if (engine == 'hdf5r') {
    cnames <- gsub(
      '/pos', '',
      gsub('cells/','',
           grep('/pos',
                grep("cells/", hdf5r::list.datasets(f), value=T),
                value=T)))
  }
  if(!is.null(cell.clusters)) {
    # count abundancies per element for each cell cluster
    cell.clusters <- as.factor(cell.clusters)
    if(!any(names(cell.clusters) %in% cnames)) {
      warning(paste("could not match any of the specified cell names. hdf5 file contains names like [",paste(cnames[1:3],collapse=' '),"... ]"))
      cat("parsing out feature counts across all cells ... ")
      if (engine == 'h5') {
        info$cluster.feature.counts <- cbind('all'=tabulate(unlist(lapply(
          cnames,
          function(n) f[paste('/cells',n,'ixs',sep='/')][] ))+1,nbins=length(info$chrm)))
      } else if (engine == 'hdf5r') {
        info$cluster.feature.counts <- cbind('all'=tabulate(unlist(lapply(
          cnames,
          function(n) f[[paste('cells',n,'ixs',sep='/')]]$read() ))+1,
          nbins=length(info$chrm)))
      } else {stop('Unknown engine. Use either hdf5r or h5.\n')}
      cat("done\n")
    } else {
      cat("parsing out info for",length(levels(cell.clusters)),"clusters: [");
      if (engine == 'h5'){
        cluster.feature.counts <- do.call(cbind,tapply(names(cell.clusters),as.factor(cell.clusters),function(ii) {
          cat(".")
          tabulate(unlist(lapply(
            ii,
            function(n) f[paste('/cells',n,'ixs',sep='/')][] ))+1,
            nbins=length(info$chrm))
        }))
      }
      if (engine == 'hdf5r') {
        cluster.feature.counts <- do.call(cbind,tapply(names(cell.clusters),as.factor(cell.clusters),function(ii) {
          cat(".")
          tabulate(unlist(lapply(
            ii,
            function(n) f[[paste('cells',n,'ixs',sep='/')]]$read() ))+1,
            nbins=length(info$chrm))
        }))
      }

      cat(". ]. done\n")
      info$cluster.feature.counts <- cluster.feature.counts;
    }
  } else {
    # combine counts on all cells
    cat("parsing out feature counts across all cells ... ")
    if (engine == 'h5'){
      info$cluster.feature.counts <- cbind('all'=tabulate(unlist(lapply(
        cnames,
        function(n) f[paste('/cells',n,'ixs',sep='/')][] ))+1,
        nbins=length(info$chrm)))
    }
    if (engine == 'hdf5r'){
      info$cluster.feature.counts <- cbind('all'=tabulate(unlist(lapply(
        cnames,
        function(n) f[[paste('cells',n,'ixs',sep='/')]]$read() ))+1,
        nbins=length(info$chrm)))
    }
    cat("done\n")
  }
  if (engine == 'h5' ) h5::h5close(f)
  if (engine == 'hdf5r') f$close_all()

  # calculate dataset-wide effective gene length and other parameters
  # attempt to get unique gene names
  #gchr <- paste(info$features_gene,info$chrm,sep=':')
  genes <- info$features_gene;

  if(!is.null(internal.priming.info)) {
    internal.priming.info$chr <- gsub("^chr","",as.character(internal.priming.info$chr))
  }

  t.get.lengthinfo <- function(fc,min.exon.count=10) {
    #fc <- rowSums(cluster.feature.counts)
    # find last expressed exon for each gene
    gdf <- do.call(rbind,mclapply(sn(unique(genes)),function(gene) {
      ii <- which(genes==gene);
      exons <- info$is_intron[ii]=='FALSE'
      # select exons with a valid minimal expression level
      valid.exons <- fc[ii[exons]]>=min.exon.count
      # is gene more on one chromosome, ignore the gene
      if(length(unique(info$chrm[ii]))>1) {
        valid.exons[valid.exons] <- FALSE;
      }

      if(any(valid.exons)) {
        # number of exons, their lengths
        valid.exons.sizes <- info$start_end[ii[exons][valid.exons],,drop=F]
        el <- flatLength(valid.exons.sizes[order(valid.exons.sizes[,1,drop=F],decreasing=FALSE),,drop=F]);
        # define effective range
        gene.range <- range(valid.exons.sizes)
        gene.size <- diff(gene.range)+1
        rv <- c(el=log10(el+1),il=log10((gene.size-el)+1),nex=nrow(valid.exons.sizes))
        if(!is.null(internal.priming.info)) {
          vi <- which(internal.priming.info$chr==info$chrm[ii[1]] & internal.priming.info$start>=gene.range[1] & internal.priming.info$end<=gene.range[2])
          if(length(vi)>0) {
            nipconc <- sum(internal.priming.info$conc[vi]);
            rv <- c(rv,c(nipconc=nipconc,nipdisc=length(vi)-nipconc))
          } else {
            rv <- c(rv,c(nipconc=0,nipdisc=0))
          }
        }
        return(rv)
      } else {
        return(NULL)
      }
    },mc.cores=n.cores,mc.preschedule=TRUE))
  }
  cat("calculating gene stats ... ")
  base.df <- t.get.lengthinfo(rowSums(info$cluster.feature.counts),min.exon.count = min.exon.count)
  cat("done\n")
  return(list(gene.df=data.frame(base.df),info=info))
}

# estimate projected delta given x'=(y-o) - gamma*x solution
# em - normalized expression matrix
# nm - normalized nascent matrix
# gamma - inferred degradation coefficients
# o - inferred offset (assumed to be zero by default)
# delta - time to project forward
t.get.projected.delta <- function(em,nm,gamma,offset=rep(0,length(gamma)),delta=0.5) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em));
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn];
  # time effect constant
  egt <- exp(-gamma*delta);
  y <- nm-offset; y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  em*egt + (1-egt)*y/gamma  - em
}

# conservative estimate based on the sensitivity to addition/subtraction of a single n count
t.get.projected.delta2 <- function(em,nm,nm.size,gamma,offset=rep(0,length(gamma)),delta=0.5,scount=1) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em));
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn];
  # time effect constant
  egt <- exp(-gamma*delta);
  y <- nm-(offset %o% nm.size);
  y1 <- y-scount; y2 <- y+scount;
  y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  y1[y1<0] <- 0; y2[y2<0] <- 0;
  d <- em*egt + (1-egt)*t(t(y)/nm.size)/gamma  - em;
  d1 <- em*egt + (1-egt)*t(t(y1)/nm.size)/gamma  - em;
  d2 <- em*egt + (1-egt)*t(t(y2)/nm.size)/gamma  - em;
  
  cd <- d;
  zi <- abs(cd)>abs(d1); cd[zi] <- d1[zi]
  zi <- abs(cd)>abs(d2); cd[zi] <- d2[zi]
  cd[sign(d1)!=sign(d2)] <- 0;
  cd
}


# estimate projected delta given log2 fold observed/expected nascent ratio 
# em - normalized expression matrix
# nm - normalized nascent matrix
# gamma - inferred degradation coefficients
# r - log2 observed/expected nascent cont ratio
# delta - time to project forward
t.get.projected.delta.from.log2ratio <- function(em,gamma,r,delta=0.5,min.val=1e-4) {
  # adjust rownames
  gn <- intersect(intersect(names(gamma),rownames(em)),rownames(r));
  em <- em[gn,]; gamma <- gamma[gn]; r <- 2^r[gn,];
  # time effect constant
  egt <- exp(-gamma*delta);
  (em+min.val)*(egt*(1-r) +r) - em
}


# determine membership of points in fragments
points.within <- function(x,fs,fe,return.list=F,return.unique=F,sorted=F,return.point.counts=F) {
  if(is.null(x) | length(x) < 1) { return(c()) };
  if(!sorted) {
    ox <- rank(x,ties.method="first");
    x <- sort(x);
  }
  
  se <- c(fs,fe);
  fi <- seq(1:length(fs));
  fi <- c(fi,-1*fi);
  
  fi <- fi[order(se)];
  se <- sort(se);
  
  storage.mode(x) <- storage.mode(fi) <- storage.mode(se) <- "integer";
  if(return.unique) { iu <- 1; } else { iu <- 0; }
  if(return.list) { il <- 1; } else { il <- 0; }
  if(return.point.counts) { rpc <- 1; } else { rpc <- 0; }
  storage.mode(iu) <- storage.mode(il) <- storage.mode(rpc) <- "integer";
  result <- points_within2(x,se,fi,il,iu,rpc)
  #result <- .Call("points_within2",x,se,fi,il,iu,rpc);
  if(!sorted & !return.point.counts) {
    result <- result[ox];
  }
  return(result);
}



balancedKNN <- function(val,k,maxl=k,return.distance.values=FALSE,n.threads=1,dist='cor') {
  if(class(dist)=="dist") { # actual distance was passed
    if(!all(labels(dist)==colnames(val))) { stop("balancedKNN(): supplied distance doesn't match the columns of val") }
    cd <- as.matrix(dist);
  }  else {
    if(dist=='cor') {
      cd <- 1-cor(val);
    } else if(dist=='euclidean') {
      cd <- as.matrix(dist(t(val)))
    } else {
      stop(paste("unknown distance",dist,"specified"))
    }
  }
  z <-  balanced_knn(cd,k,maxl,return.distance.values,n.threads);
  rownames(z) <- colnames(z) <- colnames(val);
  z
}


# fater matrix correlations wtih armadillo
##' A slightly faster way of calculating column correlation matrix
##' @param mat matrix whose columns will be correlated
##' @param nthreads number of threads to use 
##' @return correlation matrix 
##' @export
armaCor <- function(mat,nthreads=1) {
  cd <- arma_mat_cor(mat);
  rownames(cd) <- colnames(cd) <- colnames(mat);
  return(cd)
}



defaultNCores <- function() { parallel::detectCores(logical=F) }
