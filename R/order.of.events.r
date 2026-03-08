#' Leveraging bnlearn for inferring the order of complex events.
#' This function uses the bnlearn package to infer the order of complex events based on provided data.
#' It constructs a Bayesian network using the given complex events, breakpoints, and coverage data to determine their relationships.
#' @param cplx_list A list containing information about complex events, including blocks and contigs.
#' @param bpdf A data frame containing breakpoint information, with columns for contig names, start positions, and end positions.
#' @param covdf A data frame containing coverage information for the contigs.
#' @param cplx.event The index of the complex event to be analyzed.
#' @param optimization The type of optimization to perform ("hc" for hill climbing, "pc.stable" for PC-stable algorithm).
#' @importFrom reshape2 dcast
#' @importFrom dplyr mutate
#' @export
order.of.events <- function(cplx_list, bpdf, covdf, cplx.event = 1,optimization = "pc.stable") {
    library(bnlearn)
    library(reshape2)
    library(dplyr)
    i = cplx.event
    t.cplx.blk = cplx_list$blocks[i]
    t.cplx.ctg = cplx_list$contigs[i]
    t.cplx.bp = bpdf[which(bpdf[, 1] %in% unlist(t.cplx.ctg)), ]
    #rm(cplx_list)
    
    cplx.bp.colnames = unique(as.vector(paste0(t.cplx.bp[, 2], ":", t.cplx.bp[, 3])))
    cplx.bp.rownames = unique(t.cplx.bp$V1)
    
    ## generating the contig vs breakpoint matrix
    t.cplx.bp$V4 <- paste(paste0(t.cplx.bp$V2, "0kb"), paste0(t.cplx.bp$V3, "0kb"), sep = ":")
    t.cplx.bp <- t.cplx.bp[, c("V1", "V4")]
    t.cplx.bp$value = 1
    cplx.bp.mtx <- reshape2::dcast(t.cplx.bp, V1 ~ V4)
    if(nrow(covdf)>0){
    cplx.bp.mtx.cov=merge(cplx.bp.mtx,covdf, by.x="V1", by.y="mapID")

   ## Adding contig coverage to reinforce evidence for each event type
    cplx.bp.mtx.rep <- as.data.frame(lapply(cplx.bp.mtx.cov, rep, round(cplx.bp.mtx.cov$cov)))
    cplx.bp.mtx.rep=subset(cplx.bp.mtx.rep, select=-c(V1,cov))
    colnames( cplx.bp.mtx.rep)<-gsub("_",".",colnames( cplx.bp.mtx.rep))
    cplx.bp.mtx.rep[, 1:ncol(cplx.bp.mtx.rep)] <- sapply(cplx.bp.mtx.rep[, 1:ncol(cplx.bp.mtx.rep)], as.numeric) 
    } else {
  	cplx.bp.mtx.rep<-subset(cplx.bp.mtx, select=-V1)
  	cplx.bp.mtx.rep[, 1:ncol(cplx.bp.mtx.rep)] <- sapply(cplx.bp.mtx.rep[, 1:ncol(cplx.bp.mtx.rep)], as.numeric)
    } 

    ## Bayesian optimization step to get the tree object.
    if (optimization == "hc") {
      res = hc(cplx.bp.mtx.rep)
    }
    if (optimization == "pc.stable") {
      res = pc.stable(cplx.bp.mtx.rep)
      
      ## remove nodes without any ingoing and outgoing dirction
      empty_nodes <- c()
      for (nodes in 1:length(res$nodes)) {
        if (length(unlist(res$nodes[[nodes]])) == 0)
        {
          empty_nodes <- c(empty_nodes, nodes)
        }
      }
      if (length(empty_nodes) < length(res$nodes)) {
        modified_list <- res$nodes[setdiff(seq_along(res$nodes), empty_nodes)]
        res$nodes <- modified_list
      } else NULL
    }
    return(res)
  }
