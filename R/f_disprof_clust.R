#' Dissimilarity profile analysis (DISPROF) of a cluster analysis dendrogram
#'
#' This function follows Clarke et al.'s (2008) description of similarity
#' profile analysis (SIMPROF), but employs an equivalent approach based on
#' dissimilarity profile analysis (DISPROF).\cr
#' \cr
#' This method is a form of agglomerative, hierarchical' cluster analysis.
#'
#' @references Clarke, K. R., P. J. Somerfield, and R. N. Gorley. 2008.  Testing
#'   null hypotheses in exploratory community analyses: similarity profiles and
#'   biota-environmental linkage. J. Exp. Mar. Biol. Ecol. 366:56-69.
#'   (\href{https://doi.org/10.1016/j.jembe.2008.07.009}{Clarke et al., 2008})
#' @references Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English
#'   ed. Elsevier Science BV, Amsterdam.
#'
#' @param Y    Matrix of response data (rows = observations, columns = variables)
#' @param dis  Dissimilarity measure to be applied (see \code{vegan::vegdist}
#'   for additional details):
#' \itemize{
#'  \item Euclidean           (\code{'euc'})
#'  \item Bray Curtis         (\code{'bray'})
#'  \item Manhattan           (\code{'manhattan'})
#'  \item Gower's             (\code{'gower'})
#'  \item Alternative Gower's (\code{'altGower'})
#'  \item Canberra            (\code{'canberra'})
#'  \item Clark               (\code{'clark'})
#'  \item Kulczynski          (\code{'kulczynski'})
#'  \item Morisita            (\code{'morisita'})
#'  \item Horn                (\code{'horn'})
#'  \item Binomial            (\code{'binomial'})
#'  \item Cao                 (\code{'cao'})
#' }
#'
#' @param txt  Vector of row labels; generic labels are automatically generated if not provided
#' @param link Type of cluster analysis the be performed:
#' \enumerate{
#'  \item average    (UPGMA)
#'  \item centroid   (UPGMC)
#'  \item median     (WPGMC)
#'  \item single     (Shortest Distance)
#'  \item ward.D2    (Inner Squared Distance)
#'  \item mcquitty   (WPGMA)
#' }
#' @param mc Toggle progressive adjustment of p-value for multiple comparisons
#' @param iter Number of permutations for DISPROF analysis
#' @param alpha Significance level (default = 0.05)
#' @param tol Tolerance for rejecting a p-value > alpha (default = 0.005)
#'
#' @return Returns a list structure that contains the following:
#' \itemize{
#'  \item \code{grp}: vector specifying group membership for each row of Y
#'  \item \code{inc}: record of incremental change made to GRP
#'  \item \code{p}: corresponding randomized p-values
#'  \item \code{p_idx}: corresponding rows of Y used in the assessment
#'  \item \code{adj}: cell array of string indicating whether p-values were
#'  adjusted
#'  \item \code{Z}: cluster analysis tree linkages
#'  \item \code{cols}: vector indicating splits that were associated with
#'  significant p-values (=1) and those that were note (=0)
#'  \item \code{method}: method used for constructing cluster analysis tree
#'  \item \code{bin}: symmetric binary connectivity matrix defining pair-wise
#'  relationships as objects that are members of the same cluster (=0) or
#'  members of different clusters (=1)'
#' }
#'
#'
#' @export

f_disprof_clust <- function(Y,
                            dis    = 'euc',
                            txt    = NULL,
                            link   = 1,
                            mc     = FALSE,
                            iter   = 1000,
                            alpha  = 0.05,
                            tol    = 0.005){

# Checks user defined link type to pass to hclust().
if        (link == 1) {
  linkTxt  <- 'average'
  linkType <- 'UPGMA'
} else if (link == 2){
  linkTxt  <- 'centroid'
  linkType <- 'UPGMC'
} else if (link == 3){
  linkTxt  <- 'median'
  linkType <- 'WPGMC'
} else if (link == 4){
  linkTxt  <- 'single'
  linkType <- 'Shortest Distance'
} else if (link == 5){
  linkTxt  <- 'ward.D2'
  linkType <- 'Inner Squared Distance'
} else if (link == 6){
  linkTxt  <- 'mcquitty' # weighted
  linkType <- 'WPGMA'
} else {
  warning("'link' must be a value from 1 to 6'")
}

# ensure Y is in appropriate format for distance measures
Y <- as.matrix(Y)


# Set default grouping vector
if (is.null(txt) == TRUE) {
  txt <- as.character(1:(dim(Y)[1]))
}

# Check whether labels are of compatible size
if (length(txt) != dim(Y)[1]) {
  warning("'txt' does not have the same number of rows as 'Y'")
} else {
  rownames(Y) <- txt
}

# Get number of rows of input data
nR = dim(Y)[1]

################################################################################
#                            CLUSTER ANALYSIS TREE:                            #
################################################################################

# Perform agglomerative heirarchical cluster analysis.
# extracts the lower tridiagonal of the the resemblance matrix of Y, based on
# user defined method 'dis', and method 'linkTxt'
Y_dis     <- vegan::vegdist(Y, dis, na.rm = TRUE)

Z         <- hclust(Y_dis,linkTxt)

# Extract cut-off values, descending order
cv        <- rev(Z$height)
# Replaces any distances = zero with eps() in case duplicate rows occur in Y)
cv[cv==0] <- pracma::eps(1)


### Determine the cluster membership at each split
cluT <- matrix(NaN, nR, nR-1)
for (i in 1:length(cv)) {
  T         <- cutree(Z, h = cv[i])
  cluT[,i]  <- f_dummy2cat(as.matrix(vegan::vegdist(T, 'euc')) == 0)
}

### Recode so cluster numbers start with 1:
oldC <- unique(as.vector(cluT)) # old cluster number
nClu <- length(oldC)            # get no. of clusters
newC <- 1:length(oldC)          # new cluster number
clu  <- matrix(NaN, nR, nR-1)   # pre-allocate null matrix)
for (i in 1:nClu) {
  clu[cluT==oldC[i]] <- newC[i]
}
if (identical(clu[,1], rep(1,nR)) == FALSE) {
  warning("Column 1 of CLU doesn't specify a SINGLE group!")
}

# Create binary version of 'clu' so 0's indicate clusters that haven't
# changed since the previous split (col) and 1's indicate those that have; this
# is used to avoid assessing the same cluster more than once by DISPROF:
cluN <-  clu                      # make a copy
uGrp <-  unique(as.vector(cluN))  # get unique groups
nGrp <-  length(uGrp)             # get no. of groups
for (i in 1:nGrp) {
  idxCol = which(apply((abs(t(diff(t(cluN==uGrp[i])*1)))),2,sum)==0)
  idxCol = idxCol+1
  for (j in idxCol) {
    cluN[cluN[,j]==uGrp[i],j] <- 0
  }
}

# Takes logical index, turns it into binary matrix for subsequent math
cluB <- (cluN>0)*1

# Cleanup environment
rm(i, j, oldC, newC, cluT, cluN)

################################################################################
#                Test of HOMOGENEITY OF COMPOSITION via DISPROF                #
################################################################################
nC    <- dim(clu)[2]   # get no. of columns
grp   <- rep(1, nR)    # Initialize as single group
grpB  <- rep(1, nR)    # Initialize as single group
inc   <- matrix(NaN, nR, nC-1)
Pi    <- rep(NaN, nR)
p     <- rep(NaN, nR)
p_idx <- list()
cnt   <- 0             # Initialize p-value counter
cntG  <- 1             # Initial group counter
colS  <- rep(0, nC)    # Initialize significant col indicator

cat('\n========================================================
     \n              Cluster Analysis via DISPROF
     \n--------------------------------------------------------
     \nEvaluating presence of...')

for (i in 1:(nC-1)) {
  idx       <- grp > 0        # get index of rows not in terminal groups
  grp[idx]  <- clu[idx,i]     # update rows not marked 'terminal
  grpB[idx] <- cluB[idx,i+1]  # update corresponding binary vector
  inc[,i]   <- grp            # record incremental change to grouping vector

  # Get sorted list of unique groups, but flagged with grpB
  uGrp      <- sort(unique(grp*grpB))

  # Remove those marked as terminal groups (< 0) or haven't changed since the
  # previous column (= 0)
  uGrp <- uGrp[uGrp > 0]

  if (length(uGrp) != 0) {
    nGrp <- length(uGrp)          # get no. of group (there should only be 1)
    for (j in 1:nGrp) {           # Process each group:

      idxG <- which(grp==uGrp[j]) # Get index to rows of this group

      # Check if min no. of rows are present:
      if (length(idxG)<2) {
        # Mark rows as terminal group with a negative number
        grp[idxG] <- grp[idxG] * -1
      } else {
        # Check for the presence of multivariate structure in this group:
        cnt <- cnt + 1
        cat("\n", cntG+1, 'groups: ')

        # check rows for structure
        disprof      <- f_disprof(Y[idxG,], dis, iter, FALSE, FALSE, FALSE)
        Pi[cnt]      <- disprof$Pi  # Collect Pi Statistics
        p[cnt]       <- disprof$p   # Collect p-values
        p_idx[[cnt]] <- idxG        # Collect list of rows tests

        # -- Optionally adjust p-value for multiple comparisons --
        if (mc > 0) {
          p[cnt] <- p[cnt]*cnt #Progressive Bonferroni method
          if (p[cnt] > 1) {
            p[cnt] = 1
          }
        }
        cat(paste0("Pi = ", format(round(Pi[cnt],3), nsmall = 4), ", ",
            "p = ", format(round(p[cnt],3), nsmall = 4), sep = ""))

        # If not significant structure exists, no further splitting of group
        if ((p[cnt] - alpha)>tol) {
          # Mark rows as terminal group with a negative number
          grp[idxG] <- grp[idxG] * -1
        } else {
          cntG    <- cntG + 1
          colS[i] <- colS[i] + 1
        }
      }
    }
  }
}

if (mc > 0) {
  cat("\n--------------------------------------------------------",
      "\n No. of groups identified      = ", cntG,
      "\n No. of permutation iterations = ", iter,
      "\n alpha level                   = ", alpha,
      "\n dissimilarity metric          = ", dis,
      "\n p-values adjusted for multiple testing: YES",
      "\n--------------------------------------------------------\n")
} else {
  cat("\n--------------------------------------------------------",
      "\n No. of groups identified      = ", cntG,
      "\n No. of permutation iterations = ", iter,
      "\n alpha level                   = ", alpha,
      "\n dissimilarity metric          = ", dis,
      "\n p-values adjusted for multiple testing: NO",
      "\n--------------------------------------------------------\n")
}

# format for output (remove flag indicating terminal groups)
grp <- abs(grp)
inc <- abs(inc)

# Recode grouping vector as consecutive integers
grpOld  <- grp             # make a copy
incOld  <- inc
grp     <- f_recode(grp)   # recode as consecutive integers
uGrpOld <- unique(grpOld)  # get unsorted list of original groups
uGrp    <- unique(grp)     # get unsorted list of new groups
nGrp    <- length(uGrpOld) # get no. of groups

for (i in 1:nGrp) {
  inc[incOld==uGrpOld[i]] <- uGrp[i]
}

# Remove NaNs
idxNaN <- is.nan(p)
p      <- p[!idxNaN]
p_idx  <- p_idx[!idxNaN]

# References 'grp' value to build lists of member IDs
clusters         <- list()
for (i in 1:length(unique(grp))) {
  clusters[[i]] <- txt[grp == i]
}

# Wrap results into output structure
result          <- list() # Pre-allocate the list

result$grp      <- grp
result$clusters <- clusters
result$inc      <- inc
result$Pi       <- Pi
result$p        <- p
result$p_idx    <- p_idx
if (mc > 0) {
  result$adj    <- "p-values adjusted"
} else {
  result$adj    <- "p-values not adjusted"
}
result$Z        <- Z #list of cluster linkages needed to draw dendrograms
result$colS     <- colS
result$linkType <- linkType

    ### MS COMMENTS (2022-01-25): Need to reconcile this output - in Dr. Jones
    ### version, when 2 outputs are defined (i.e., [result,bin]), then it
    ### performs an additional step). Maybe just force this as part of the
    ### output structure and call it a day?

      # Optionally return binary connectivity matrix
      #   0 = objects in same cluster
      #   1 = objects in different clusters
      # if(nargout > 1){
      # bin <- ~(f_dis(grp, 'euc') ==0)
      #}
  return(result)
}
