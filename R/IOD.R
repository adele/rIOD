#' @export IOD
IOD <- function(suffStat, alpha=0.05, method = "standard", procedure = "original", verbose=FALSE) {

  initSkeletonOutput <- initialSkeleton(suffStat, alpha, procedure=procedure, verbose=verbose)
  G <- initSkeletonOutput$G
  # renderAG(G)

  IP <- initSkeletonOutput$IP
  sepsetList <- initSkeletonOutput$sepsetList
  # lapply(sepsetList, formatSepset)

  Gi_list <- initSkeletonOutput$listGi
  # lapply(Gi_list, renderAG)

  nCK1 <- sum(unlist(initSkeletonOutput$nCK1List), na.rm = TRUE)
  nNCK <- sum(unlist(initSkeletonOutput$nNCKList), na.rm = TRUE)

  n_datasets <- length(suffStat$citestResultsList)

  possImmfromTriplets <- initSkeletonOutput$possImmfromTriplets

  # Algorithm 3:

  p <- length(colnames(G))

  G_copy <- G
  G_copy[which(G_copy == 3)] <- 1

  #all((amat != 0) == (t(amat != 0))) ist nicht TRUE
  possSepList <- setOfPossSep(G_copy,p)

  existingEdges <- adjPairsOneOccurrence(G)

  RemEdges <- getRemEdges(existingEdges,G, possSepList,n_datasets,suffStat)


  power_RemEdges <- powerSet(unique(RemEdges))
  # one_edge_list <- which(lapply(power_RemEdges, length) == 1)

  index_possImmList <- 1

  #G_PAG_List <- list()
  #for (E in power_RemEdges) {
  G_PAG_List <- foreach (E = power_RemEdges, .verbose=verbose) %dofuture% {
    H <- induceSubgraph(G,E)
    labelsG <- colnames(G)
    PossImm <- getPossImm(H, n_datasets, suffStat, sepsetList, labelsG)

    if (procedure == "orderedtriplets") {
      savetails <- H
      savetails[which(savetails != 3)] <- 0

      H[which(savetails == 3)] <- 1
    }

    listAllHt <- list()
    if (length(PossImm) > 0) {
      power_possImm <- all_combinations(PossImm)
      for (t in power_possImm) {
        H_t <- H
        # orient Colliders
        for (tau in t) {
          if (!is.null(tau)) {
            H_t[tau[1], tau[2]] <- 2
            H_t[tau[3], tau[2]] <- 2
          }
        }
        listAllHt[[length(listAllHt)+1]] <- H_t
      }
    } else {
      listAllHt[[length(listAllHt)+1]] <- H
    }


    G_PAG <- unique(applyRulesOnHt(unique(listAllHt)))
    temp_G_PAG <- G_PAG
    # TODO: below, just change circles to tails, not arrowhead to tails..
    # if there is an arrowhead on a place that a definite tail is expected,
    # we should remove the PAG inside of this loop
    # -- it is fine to remove these PAGs before
    # getting the G_PAG_List_before.


    if(procedure == "orderedtriplets"){
      #include Tails
      indices <- which(savetails == 3, arr.ind = TRUE)
      # Zeilennamen der Werte, die 3 sind
      row_names <- rownames(savetails)[indices[, 1]]

      # Spaltennamen der Werte, die 3 sind
      col_names <- colnames(savetails)[indices[, 2]]

      if (length(row_names) > 0) {
        temp_G_PAG <- list()
        for (i in 1:length(G_PAG)) {
          exclude_Pag <- FALSE
          #% should not be indices but names
          cur <- G_PAG[[i]]

          for (j in 1:length(row_names)) {
            # do not want to force the arrowheads
            if (cur[row_names[[j]], col_names[[j]]] != 2) {
              cur[row_names[[j]], col_names[[j]]] <- 3
            } else{
              exclude_Pag <- TRUE
            }
            # idk what I did here, maybe a left over from something that passed away
            #if(i == length(G_PAG) && j == length(row_names)){
            #exclude_Pag <- TRUE
            #}
          }

          if (!exclude_Pag) {
            temp_G_PAG[[length(temp_G_PAG)+1]] <- cur
          }
        }
      }
    }
    G_PAG <- temp_G_PAG
    # For each possible G in the power set of graphs you are creating, make sure
    # to update sepset accordingly.
    return(G_PAG)
  }

  # This is counting the number of PAGs before using the violation checks,
  # which takes more time in the computational sense.
  if (length(unique(unlist(G_PAG_List))) != 0) {
    #The if excludes list() and list(list(), list(), list()), list(list(list(), list()), list(list(), list())) ...
    G_PAG_List_before <- unique(unlist(G_PAG_List, recursive=F))
    len_before <- length(G_PAG_List_before)
  } else {
    G_PAG_List_before <- G_PAG_List <- list()
    len_before <- 0
  }

  if (len_before > 0) {
    violation_List <- validatePossPags(G_PAG_List_before, sepsetList, suffStat, IP, method, Gi_list)
    G_PAG_List <- G_PAG_List_before[!violation_List]
  }
  return(list(G_PAG_List=G_PAG_List, Gi_PAG_list=Gi_list,
              G_PAG_List_before=G_PAG_List_before, len_before=len_before,
              nCK1=nCK1, nNCK=nNCK))
}
