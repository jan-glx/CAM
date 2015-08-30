CAM <-
function(X, scoreName = "SEMGAM", 
                          parsScore = list(numBasisFcts=10), 
                          numCores = 1, 
                          maxNumParents = min(dim(X)[2] - 1, round(dim(X)[1]/20)),
                          output = FALSE,
                          fixedOrders = matrix(ncol=2,nrow=0),
                          orderFixationMethod = "force_edge",
                          variableSel = FALSE, 
                          variableSelMethod = selGamBoost, 
                          variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),
                          pruning = FALSE, 
                          pruneMethod = selGam, 
                          pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts=10),
                          intervData = FALSE, 
                          intervMat = NA) 
{
     # check whether fixedOrders input is correct
    possibleOrderFixationMethods <- c("force_edge", "emulate_edge", "free")
    if (!(orderFixationMethod %in% possibleOrderFixationMethods)) 
      stop("orderFixationMethod must be either ",paste0(possibleOrderFixationMethods, collapse = " or "),".")
    if (class(fixedOrders) == "numeric" && length(fixedOrders) == 2){
      fixedOrders <- matrix(as.integer(fixedOrders), nrow=1)
    } else if (is.null(fixedOrders)) {
      fixedOrders <- matrix(ncol=2,nrow=0)
    } else if (!is.matrix(fixedOrders) || ncol(fixedOrders) != 2){
      stop("fixedOrders must be a [nrow X 2] matrix where pi(fixedOrders[i,1]) < pi(fixedOrders[i,2]), and pi is a causal ordering of the nodes")
    }
  
    if(output)
    {
        cat("number of cores:", numCores, "\n")
    }
    # We record the time consumption. They are shown if output == TRUE
    timeCycle <- 0
    timeUpdate <- 0
    timeScoreMat <- 0
    timeSel <- 0
    timePrune <- 0
    timeMax <- 0
    
    # we record how the score develops 
    scoreVec <- c()
    # and which edges are added
    edgeList <- c()
    
    # this counter is only used if output = TRUE
    counterUpdate <- 0
    p <- dim(X)[2]
    
    
    ####
    # STEP 1: variable selection
    ####
    # A matrix selMat is constructed. Entry (i,j) being one means that i is a possible parent of j.
    if(variableSel)
    {
            ptm <- proc.time()[3]
            if(intervData)
            {
                X2 <- X[rowSums(intervMat) == 0,]
                if(output)
                    cat("The preliminary neighbourhood selection is done with the observational data only.\n")
            } else
            {
                X2 <- X
            }
            if(numCores == 1)
            {
                selMat <- mapply(variableSelMethod,MoreArgs = list(X = X2, pars = variableSelMethodPars, output = output),1:p)
            } else
            {
                selMat <- mcmapply(variableSelMethod,MoreArgs = list(X = X2, pars = variableSelMethodPars, output = output),1:p, mc.cores = numCores)
            }
            # The next line includes j as a possible parent of i if i is considered a possible parent of j
            # selMat <- selMat | t(selMat)
            cou <- 0
            for(jk in 1:p)
            {
                cou <- cou + 2^{sum(selMat[,jk])}
            }
            if(output)
            {
                cat("Instead of p2^(p-1) -Sillander- ",p*2^(p-1) ," we have ", cou, "\n")
                cat("Greedy, on the other hand, is computing ",sum(selMat) ," entries. \n")
            }
            timeSel <- timeSel + proc.time()[3] - ptm
    } else
    {
        selMat <- matrix(TRUE, p,p)
    }
    if(variableSel & output)
    {
        if(output)
        {
            if(p<30)
            {
                cat("This is the matrix of possible parents after the first step.\n")
                show(selMat)
            }
            cat("Object size of selmat: ", object.size(selMat), "\n")
        }
    }
    
    
    ####
    # STEP 2: Include Edges
    ####
    # compute score matrix 
    ptm <- proc.time()[3]
    computeScoreMatTmp <- computeScoreMat(X, scoreName=scoreName, numParents = 1, numCores = numCores, output = output, 
                                  selMat = selMat, parsScore = parsScore, intervMat = intervMat, intervData = intervData)
    timeScoreMat <- timeScoreMat + proc.time()[3] - ptm
    if(output)
    {
        cat("Object size of computeScoreMatTmp: ", object.size(computeScoreMatTmp), "\n" )
    }
    # We need the pathMatrix (entry (i,j) being one means that there is a directed path from i to j) in order to keep track of possible cycles.
    pathMatrix <- matrix(FALSE,p,p)
    diag(pathMatrix) <- TRUE
    if (orderFixationMethod == "emulate_edge"){
      computeScoreMatTmp$scoreMat[fixedOrders[,2:1,drop=F]] <- -Inf #prevents adding of a edge from j to i 
      pathMatrix[fixedOrders] <- TRUE #will prevent adding of edges from DE(j) to {AC(i),i} and {DE(j),j} to {AC(i}
      # where (i,j) is a row of fixedOrders
    }
        
    Adj <- sparseMatrix(i=integer(),j=integer(),dims=c(p,p))
    scoreNodes <- computeScoreMatTmp$scoreEmptyNodes
  
    fixedOrdersAdded <- 0L
      
    # Greedily adding edges
    while(sum(computeScoreMatTmp$scoreMat!=-Inf) > 0)
    {
        # Find the best edge
        ptm <- proc.time()[3]
        
        if (orderFixationMethod == "force_edge" && fixedOrdersAdded < nrow(fixedOrders)){
          fixedOrdersAdded = fixedOrdersAdded + 1L
          ix_max <- fixedOrders[fixedOrdersAdded, , drop=F]
        } else {
          ix_max <- arrayInd(which.max(computeScoreMatTmp$scoreMat), dim(computeScoreMatTmp$scoreMat))
        }
        
        timeMax <- timeMax + proc.time()[3] - ptm
        Adj[ix_max] <- TRUE
        scoreNodes[ix_max[2]] <- scoreNodes[ix_max[2]] + computeScoreMatTmp$scoreMat[ix_max]
        if(output)
        {
            cat("\n Included edge (from, to) ", ix_max, "\n")
        }
        
        # Do not include the same edge twice.
        computeScoreMatTmp$scoreMat[ix_max] <- -Inf
        
        # Avoid cycles
        ptm <- proc.time()[3]
        DescOfAndNewChild <- which(pathMatrix[ix_max[2],])
        AncOfAndNewParent <- which(pathMatrix[,ix_max[1]])
        pathMatrix[AncOfAndNewParent,DescOfAndNewChild] <- TRUE
        computeScoreMatTmp$scoreMat[DescOfAndNewChild,AncOfAndNewParent] <- -Inf 
        timeCycle <- timeCycle + proc.time()[3] - ptm
        
        # Record the score of the current graph
        scoreVec <- c(scoreVec, sum(scoreNodes))
        # Record which edge has been added
        edgeList <- rbind(edgeList, ix_max, deparse.level=0)
        
        # Update the score of column j
        ptm <- proc.time()[3]
        computeScoreMatTmp$scoreMat <- updateScoreMat(computeScoreMatTmp$scoreMat, X, scoreName = scoreName, ix_max[1], ix_max[2],
                                                      scoreNodes, Adj, numCores=numCores, output = output, maxNumParents = maxNumParents, 
                                                      parsScore = parsScore, intervMat = intervMat, intervData = intervData)
        timeUpdate <- timeUpdate + proc.time()[3] - ptm
        
        counterUpdate <- counterUpdate + 1
    }
    
    
    ####
    # STEP 3: Prune the DAG
    ####
    if(pruning)
    {
        if(intervData)
        {
            X2 <- X[rowSums(intervMat) == 0,]
            cat("The preliminary neighbourhood selection is done with the observational data only.\n")
        } else
        {
            X2 <- X
        }
        if(output)
        {
            cat("\n Performing pruning ... \n ")
        }
        ptm <- proc.time()[3]
        Adj <- pruning(X=X2,G=Adj,pruneMethod = pruneMethod, pruneMethodPars = pruneMethodPars, output=output)      
        timePrune <- timePrune + proc.time()[3] - ptm          
    }
    
    
    
    ####
    # Output and return
    ####
    timeTotal <- timeSel + timeScoreMat + timeCycle + timeUpdate + timeMax + timePrune
    if(output)
    {
        cat("amount of time for variable selection:",timeSel,"\n")
        cat("amount of time computing the initial scoreMat:",timeScoreMat,"\n")
        cat("amount of time checking for cycles:",timeCycle,"\n")
        cat("amount of time computing updates for the scoreMat:",timeUpdate,", doing",counterUpdate,"updates.\n")
        cat("amount of time for pruning:",timePrune,"\n")
        cat("amount of time for finding maximum:",timeMax,"\n")
        cat("amount of time in total:",timeTotal,"\n")
    }
    
    result <- list(Adj = Adj, Score = sum(scoreNodes), timesVec = c(timeSel, timeScoreMat, timeCycle, timeUpdate, timePrune, timeMax, timeTotal), scoreVec = scoreVec, edgeList = edgeList)
    return(result)  
}
