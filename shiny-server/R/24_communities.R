#### Community detection ####
# Grouping the genes a correlation network draws, and ringing each group on the
# diagram.
#
# Two things are being chosen separately, and the Network tab asks for both:
#
#   the similarity   decides which pairs are joined by an edge, and is what the
#                    adjacency list reports. Spearman, Pearson, Kendall or cosine.
#   the distance     decides which genes are near each other for the purpose of
#                    forming communities. Euclidean, manhattan or jaccard, over
#                    the gene x screen score matrix, exactly as the dendrograms on
#                    the Heatmap and Adjacency tabs use it.
#
# Keeping them apart is what lets a network built from Spearman edges be grouped
# by, say, jaccard distance. Not every method can take a distance directly:
# k-means and the Gaussian mixture need coordinates, so a non-euclidean distance
# reaches them through a classical MDS embedding, and the factorisation needs a
# non-negative matrix, so it is given the affinity 1 - d/max(d). Each is said
# plainly in the method's own comment below.
#
# Louvain and Leiden are the exception: they are graph methods, and the graph they
# are given is the network on screen, so the similarity alone decides the
# communities and no distance is involved. They also find their own number of
# communities, so they take a resolution where the others take a k. The Network
# tab shows whichever of the two controls applies.
#
# Hard methods (Agnes, Clara, k-means, Louvain, Leiden) put each gene in one
# community. Soft methods (NMF, GMM) give each gene a membership in every
# community, so a gene can sit in two and the ellipses overlap. Both shapes come
# back as a membership matrix -- an indicator for the hard methods -- so
# everything downstream has one case to handle.

#### Inputs the methods need ####

#' Coordinates for the methods that cannot take a distance.
#'
#' Euclidean distance is already the geometry of the score matrix, so the matrix
#' is its own embedding. Manhattan and jaccard are not, so they are embedded by
#' classical multidimensional scaling, which places the genes so that ordinary
#' distances between the placements approximate the chosen ones.
#'
#' @param maxDims Cap on the embedding's dimensionality. The mixture model needs
#'   one: a diagonal Gaussian per community over several hundred screens, fitted
#'   to a few hundred genes, has more parameters than data and will chase a
#'   singularity rather than a cluster.
communityCoordinates <- function(m, metric, maxDims = NULL) {
  n <- nrow(m)
  dims <- ncol(m)
  if (!is.null(maxDims)) dims <- min(dims, maxDims)
  dims <- max(2L, min(dims, n - 1L))

  if (identical(metric, "euclidean") && is.null(maxDims)) return(as.matrix(m))

  d <- craveDist(m, metric)
  coords <- tryCatch(suppressWarnings(stats::cmdscale(d, k = dims)),
                     error = function(e) NULL)
  # cmdscale drops dimensions whose eigenvalue is not positive, which a
  # non-euclidean distance can easily produce; two is the fewest that still
  # describes a plane, and below that the embedding has told us the distances
  # carry no structure worth clustering.
  if (is.null(coords) || !is.matrix(coords) || ncol(coords) < 2) {
    coords <- tryCatch(suppressWarnings(stats::cmdscale(d, k = 2)),
                       error = function(e) NULL)
  }
  if (is.null(coords) || !is.matrix(coords) || ncol(coords) < 2) return(NULL)
  rownames(coords) <- rownames(m)
  coords
}

#' A non-negative affinity matrix, for the factorisation.
#'
#' Near genes get a value close to 1 and far genes close to 0, which is what the
#' factorisation is looking for; the distance chosen in the UI is therefore what
#' the communities are built from here too.
communityAffinity <- function(m, metric) {
  d <- as.matrix(craveDist(m, metric))
  mx <- max(d, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(NULL)
  a <- 1 - d / mx
  a[!is.finite(a)] <- 0
  diag(a) <- 1
  a
}

#' The drawn network, as a weighted undirected graph.
#'
#' For Louvain and Leiden, which partition a graph rather than a set of points.
#' The adjacency list holds both triangles, so each pair is reduced to one edge
#' first -- modularity is unchanged by scaling every weight alike, but a graph
#' with every edge twice is not the graph the tab draws.
#'
#' Weighted by the absolute similarity, as the Fruchterman-Reingold layout on the
#' same tab already is: a strongly anti-correlated pair is as firmly tied as a
#' strongly correlated one, which for screen scores is usually what is wanted, and
#' is in any case what the picture has always shown.
#'
#' @param genes The vertex set, in the order the membership matrix will use. Given
#'   explicitly so that the rows that come back line up with the score matrix even
#'   if a gene has no edge left after the cutoff.
communityGraph <- function(edges, genes) {
  pairs <- edges %>%
    filter(from != to, !is.na(corr)) %>%
    mutate(.a = pmin(from, to), .b = pmax(from, to)) %>%
    distinct(.a, .b, .keep_all = TRUE) %>%
    transmute(from = .a, to = .b, weight = abs(corr))
  if (nrow(pairs) == 0) return(NULL)
  igraph::graph_from_data_frame(pairs, directed = FALSE, vertices = genes)
}

#### The two soft fits ####

#' Symmetric non-negative matrix factorisation of an affinity matrix.
#'
#' Finds H >= 0 with A approximately H H', by the multiplicative update of Ding,
#' He and Simon. Row i of H says how much gene i belongs to each community, and
#' because nothing constrains a row to one column, a gene can belong to several --
#' which is the point of using it here.
#'
#' The update is damped (the 0.5 + 0.5 * ratio form) rather than applied outright.
#' The undamped version is the one usually written down and it oscillates on
#' matrices like these, where every entry is positive and the diagonal is 1.
symmetricNMF <- function(A, k, maxIter = CRAVE_NMF_MAX_ITER,
                         tol = CRAVE_COMMUNITY_TOL, seed = CRAVE_COMMUNITY_SEED) {
  n <- nrow(A)
  H <- withSeed(seed, matrix(stats::runif(n * k, 0.2, 1), nrow = n, ncol = k))
  previous <- Inf
  iterations <- 0L
  for (it in seq_len(maxIter)) {
    iterations <- it
    numerator   <- A %*% H
    denominator <- H %*% crossprod(H)
    H <- H * (0.5 + 0.5 * (numerator / pmax(denominator, 1e-10)))
    H[!is.finite(H)] <- 0
    error <- sum((A - tcrossprod(H))^2)
    if (is.finite(previous) && abs(previous - error) < tol * max(1, previous)) {
      previous <- error
      break
    }
    previous <- error
  }
  list(H = H, error = previous, iterations = iterations)
}

#' Fit a Gaussian mixture by expectation maximisation, diagonal covariances.
#'
#' Diagonal rather than full: a full covariance per community needs d(d+1)/2
#' parameters, and the embedding handed to it has more dimensions than a gene set
#' has genes to spare. Variances are floored at a thousandth of each dimension's
#' overall variance, which is what stops a community of two identical genes
#' collapsing onto them with infinite likelihood -- the classic way an
#' unconstrained EM run ends.
#'
#' Started from k-means rather than at random, so a fit is reproducible and does
#' not need restarts to avoid the worst local optima.
gaussianMixture <- function(X, k, maxIter = CRAVE_GMM_MAX_ITER,
                            tol = CRAVE_COMMUNITY_TOL, seed = CRAVE_COMMUNITY_SEED) {
  n <- nrow(X)
  start <- withSeed(seed, stats::kmeans(X, centers = k, nstart = 5, iter.max = 50))
  mu <- start$centers
  proportion <- as.numeric(table(factor(start$cluster, levels = seq_len(k))))
  proportion <- pmax(proportion, 1) / sum(pmax(proportion, 1))

  overall <- apply(X, 2, stats::var)
  overall[!is.finite(overall) | overall <= 0] <- 1
  floorVar <- 1e-3 * overall
  variance <- matrix(rep(overall, each = k), nrow = k)

  logLik <- -Inf
  responsibility <- matrix(1 / k, nrow = n, ncol = k)
  iterations <- 0L

  for (it in seq_len(maxIter)) {
    iterations <- it

    # E step, in logs throughout: the densities underflow to zero in any
    # dimensionality worth using.
    logDensity <- matrix(0, nrow = n, ncol = k)
    for (j in seq_len(k)) {
      v <- pmax(variance[j, ], floorVar)
      deviation <- sweep(X, 2, mu[j, ], "-")^2
      quadratic <- rowSums(sweep(deviation, 2, v, "/"))
      logDensity[, j] <- log(proportion[j]) - 0.5 * sum(log(2 * pi * v)) -
        0.5 * quadratic
    }
    rowMax <- apply(logDensity, 1, max)
    rowSum <- rowMax + log(rowSums(exp(logDensity - rowMax)))
    responsibility <- exp(logDensity - rowSum)

    newLogLik <- sum(rowSum)
    if (is.finite(newLogLik) && is.finite(logLik) &&
        abs(newLogLik - logLik) < tol * max(1, abs(newLogLik))) {
      logLik <- newLogLik
      break
    }
    logLik <- newLogLik

    # M step
    weight <- pmax(colSums(responsibility), 1e-10)
    proportion <- weight / n
    mu <- (t(responsibility) %*% X) / weight
    for (j in seq_len(k)) {
      deviation <- sweep(X, 2, mu[j, ], "-")^2
      variance[j, ] <- pmax(colSums(responsibility[, j] * deviation) / weight[j],
                            floorVar)
    }
  }

  list(responsibility = responsibility, logLik = logLik, iterations = iterations)
}

#### The dispatcher ####

#' Group genes into communities.
#'
#' @param m Gene x screen score matrix, rownames the gene symbols.
#' @param method One of CRAVE_COMMUNITY_METHODS.
#' @param k Number of communities, for every method that takes one. Ignored by
#'   Louvain and Leiden, which decide for themselves how many there are.
#' @param metric One of CRAVE_DISTANCES. Ignored by Louvain and Leiden, which
#'   partition the network rather than a distance between gene profiles.
#' @param edges The adjacency list the Network tab drew, columns from/to/corr.
#'   Needed by Louvain and Leiden and unused by everything else.
#' @param resolution Granularity for Louvain and Leiden; above 1 gives more and
#'   smaller communities.
#' @return list(ok, message, method, probabilistic, graph, k, resolution, metric,
#'   genes, assignment, probability, membership). `membership` is genes x k and
#'   sums to 1 along each row; for the hard methods it is an indicator, so callers
#'   have one shape to handle. `assignment` is the strongest community and
#'   `probability` how much of the gene sits in it -- 1 throughout for a hard
#'   method. `k` is what was asked for, except under Louvain and Leiden where it
#'   is what was found.
detectCommunities <- function(m, method = CRAVE_COMMUNITY_DEFAULT,
                              k = CRAVE_CLARA_K_DEFAULT,
                              metric = CRAVE_DISTANCE_DEFAULT,
                              edges = NULL,
                              resolution = CRAVE_COMMUNITY_RESOLUTION_DEFAULT) {
  off <- function() list(ok = FALSE, off = TRUE, message = "")
  fail <- function(msg) list(ok = FALSE, off = FALSE, message = msg)

  if (is.null(method) || identical(method, CRAVE_COMMUNITY_OFF)) return(off())
  if (!method %in% CRAVE_COMMUNITY_METHODS) {
    return(fail(paste0("Unrecognised community detection method '", method, "'.")))
  }
  graphMethod <- method %in% CRAVE_COMMUNITY_GRAPH
  if (graphMethod) {
    if (!needPkg("igraph")) {
      return(fail(paste0("The igraph package is not installed, so ", method,
                         " communities cannot be detected.")))
    }
    if (identical(method, "Leiden") &&
        !exists("cluster_leiden", where = asNamespace("igraph"))) {
      return(fail(paste0("This installation of igraph (", packageVersion("igraph"),
                         ") is too old for Leiden. Upgrade igraph, or choose Louvain.")))
    }
  } else if (!needPkg("cluster")) {
    return(fail("The cluster package is not installed, so communities cannot be detected."))
  }

  m <- as.matrix(m)
  m[is.na(m)] <- 0
  genes <- rownames(m)
  n <- nrow(m)

  if (n < 3) {
    return(fail("At least three genes are needed to form communities."))
  }
  if (n > CRAVE_COMMUNITY_MAX_GENES) {
    return(fail(paste0(
      "Too many genes in the network for community detection (", n, ", maximum ",
      CRAVE_COMMUNITY_MAX_GENES, "). Raise the similarity cutoff, or name fewer ",
      "genes of interest.")))
  }

  # One of the two is asked for, never both: the graph methods find their own k,
  # and nothing else has a resolution to set.
  kk <- NA_integer_
  res <- NA_real_
  if (graphMethod) {
    metric <- NA_character_
    res <- asNumericOr(resolution, NA_real_)
    if (!is.finite(res) || res < CRAVE_COMMUNITY_RESOLUTION_MIN ||
        res > CRAVE_COMMUNITY_RESOLUTION_MAX) {
      return(fail(paste0("The resolution must be between ",
                         CRAVE_COMMUNITY_RESOLUTION_MIN, " and ",
                         CRAVE_COMMUNITY_RESOLUTION_MAX, ".")))
    }
    if (is.null(edges) || nrow(edges) == 0) {
      return(fail(paste0(method, " partitions the network itself, and this one has ",
                         "no edges to partition.")))
    }
  } else {
    kk <- suppressWarnings(as.integer(asNumericOr(k, NA_real_)))
    if (is.na(kk) || kk < 2) {
      return(fail("The number of communities (k) must be 2 or more."))
    }
    if (kk >= n) {
      return(fail(paste0(
        "k must be fewer than the ", n, " genes in the network; it is ", kk, ".")))
    }
  }

  probabilistic <- method %in% CRAVE_COMMUNITY_PROBABILISTIC
  log_info(paste0("Community detection: ", method, ", ",
                  if (graphMethod) paste0("resolution ", res)
                  else paste0("k = ", kk, ", ", metric, " distance"),
                  ", ", n, " genes."))

  assignment <- NULL
  membership <- NULL

  result <- tryCatch({
    if (startsWith(method, "Agnes")) {
      # The hierarchy is built exactly as the dendrograms are, then cut into k.
      agnesName <- agnesMethod(method)
      args <- list(x = craveDist(m, metric), method = agnesName)
      if (identical(agnesName, "flexible")) args$par.method <- CRAVE_AGNES_FLEXIBLE_PAR
      ag <- do.call(cluster::agnes, args)
      assignment <- stats::cutree(stats::as.hclust(ag), k = kk)

    } else if (identical(method, "Clara")) {
      # Given the metric by name: clara implements all three itself, jaccard
      # included, and cannot be handed a dissimilarity matrix.
      cl <- cluster::clara(m, k = kk, metric = metric, stand = FALSE,
                           samples = 50, pamLike = TRUE)
      assignment <- as.integer(cl$clustering)

    } else if (identical(method, "k-means")) {
      coords <- communityCoordinates(m, metric)
      if (is.null(coords)) stop("the distances could not be embedded", call. = FALSE)
      km <- withSeed(CRAVE_COMMUNITY_SEED,
                     stats::kmeans(coords, centers = kk, nstart = 10, iter.max = 100))
      assignment <- as.integer(km$cluster)

    } else if (identical(method, "Gaussian mixture model")) {
      coords <- communityCoordinates(m, metric, maxDims = CRAVE_GMM_MAX_DIMS)
      if (is.null(coords)) stop("the distances could not be embedded", call. = FALSE)
      fit <- gaussianMixture(coords, kk)
      membership <- fit$responsibility

    } else if (graphMethod) {
      # Modularity over the drawn edges: the only two methods here that are given
      # the network rather than the gene x screen matrix, and the only two that
      # are not told how many communities to look for.
      g <- communityGraph(edges, genes)
      if (is.null(g)) stop("the network has no edges to partition", call. = FALSE)
      weights <- igraph::E(g)$weight
      # Seeded: both sweep the vertices in a random order, so two runs of the same
      # selection would otherwise ring the genes differently.
      found <- withSeed(CRAVE_COMMUNITY_SEED, if (identical(method, "Louvain")) {
        igraph::cluster_louvain(g, weights = weights, resolution = res)
      } else {
        igraph::cluster_leiden(g, objective_function = CRAVE_LEIDEN_OBJECTIVE,
                               weights = weights, resolution = res,
                               n_iterations = CRAVE_LEIDEN_ITERATIONS)
      })
      # Renumbered to 1..k with none missing: igraph is free to leave gaps, and
      # everything below indexes the membership matrix by this number.
      assignment <- as.integer(factor(as.integer(igraph::membership(found))))
      kk <- max(assignment)

    } else if (identical(method, "Non-negative matrix factorisation")) {
      affinity <- communityAffinity(m, metric)
      if (is.null(affinity)) stop("the affinity matrix could not be formed", call. = FALSE)
      fit <- symmetricNMF(affinity, kk)
      rowTotal <- rowSums(fit$H)
      membership <- fit$H / ifelse(rowTotal > 0, rowTotal, 1)
      # A gene the factorisation gave nothing to belongs nowhere in particular;
      # say so as an even split rather than silently putting it in community 1.
      membership[rowTotal <= 0, ] <- 1 / kk

    } else {
      stop(paste0("no implementation for '", method, "'"), call. = FALSE)
    }
    TRUE
  }, error = function(e) {
    log_warn(paste0("Community detection failed: ", conditionMessage(e)))
    conditionMessage(e)
  })

  if (!isTRUE(result)) {
    return(fail(paste0("Communities could not be detected: ", result, ".")))
  }

  # Every gene on its own is a valid partition and a useless picture: it is what a
  # resolution set too high for the network does, and on a network whose edges are
  # all of similar weight it arrives abruptly. Say that, rather than ring each gene
  # separately and leave the reader to work out why.
  if (graphMethod && kk >= n) {
    return(fail(paste0(
      method, " put every one of the ", n, " genes in a community of its own at ",
      "resolution ", res, ". Lower the resolution, or lower the similarity cutoff ",
      "so that the network has more edges to group by.")))
  }

  # One shape out, whichever family produced it.
  if (is.null(membership)) {
    membership <- matrix(0, nrow = n, ncol = kk)
    membership[cbind(seq_len(n), assignment)] <- 1
  }
  membership[!is.finite(membership)] <- 0
  rownames(membership) <- genes
  colnames(membership) <- paste0("Community ", seq_len(kk))

  assignment  <- max.col(membership, ties.method = "first")
  probability <- membership[cbind(seq_len(n), assignment)]
  names(assignment)  <- genes
  names(probability) <- genes

  # An empty community is not an error -- EM and k-means both abandon one
  # sometimes -- but it is worth saying, because the user asked for k of them.
  used <- length(unique(assignment))
  if (!graphMethod && used < kk) {
    log_info(paste0("Community detection returned ", used, " non-empty communities of ",
                    kk, " asked for."))
  }

  if (graphMethod) {
    log_info(paste0(method, " found ", kk, " communities at resolution ", res, "."))
  }

  list(
    ok = TRUE, off = FALSE, message = "",
    method = method, probabilistic = probabilistic, graph = graphMethod,
    k = kk, resolution = res, metric = metric,
    genes = genes, assignment = assignment, probability = probability,
    membership = membership
  )
}

#' The line above a community diagram.
#'
#' Says which measures produced the picture, because that is what a reader needs
#' to compare it with the next one. The graph methods report what they found and
#' at what resolution; the rest report what they were asked for and over which
#' distance.
communityTitle <- function(result, similarityName, nGenes) {
  if (isTRUE(result$graph)) {
    count  <- paste0(result$k, if (result$k == 1) " community" else " communities")
    detail <- paste0("resolution ", result$resolution)
  } else {
    count  <- paste0("k = ", result$k)
    detail <- paste0(result$metric, " distance")
  }
  paste0(result$method, ", ", count, ", ", nGenes, " genes<br><sup>",
         similarityName, " edges, ", detail, "</sup>")
}

#### Drawing them ####

#' A confidence ellipse around a set of points.
#'
#' The usual construction: the covariance's eigenvectors give the axes and its
#' eigenvalues their lengths, scaled to the chi-square quantile for two degrees of
#' freedom. Communities of one or two genes have no usable covariance, and
#' collinear ones have a degenerate second axis, so every axis is floored at a pad
#' derived from the size of the layout -- without which those communities would be
#' drawn as a point or a hairline and look like a bug.
ellipsePath <- function(x, y, level = CRAVE_COMMUNITY_ELLIPSE_LEVEL,
                        pad = 0.05, points = 72L) {
  if (length(x) == 0) return(NULL)
  centre <- c(mean(x), mean(y))
  theta <- seq(0, 2 * pi, length.out = points)
  circle <- cbind(cos(theta), sin(theta))

  if (length(x) < 3) {
    spread <- if (length(x) == 1) 0 else
      0.5 * sqrt(diff(range(x))^2 + diff(range(y))^2)
    radius <- max(pad, spread + pad)
    return(data.frame(x = centre[1] + radius * circle[, 1],
                      y = centre[2] + radius * circle[, 2]))
  }

  covariance <- stats::cov(cbind(x, y))
  if (any(!is.finite(covariance))) return(NULL)
  decomposition <- eigen(covariance, symmetric = TRUE)
  radii <- sqrt(pmax(decomposition$values, 0)) *
    sqrt(stats::qchisq(level, df = 2))
  radii <- pmax(radii, pad)
  pts <- circle %*% diag(radii, nrow = 2) %*% t(decomposition$vectors)
  data.frame(x = centre[1] + pts[, 1], y = centre[2] + pts[, 2])
}

#' One ellipse per community, over the nodes that belong to it.
#'
#' @param threshold How much of a gene's membership must be in a community for it
#'   to be enclosed. The hard methods store an indicator, so any threshold at or
#'   below 1 encloses exactly the assigned genes and the parameter quietly does
#'   nothing, which is why the control is only shown for the soft ones.
communityEllipses <- function(layout, result, threshold = CRAVE_COMMUNITY_MEMBERSHIP_DEFAULT,
                              level = CRAVE_COMMUNITY_ELLIPSE_LEVEL) {
  membership <- result$membership[layout$Gene, , drop = FALSE]
  extent <- max(diff(range(layout$x)), diff(range(layout$y)), na.rm = TRUE)
  pad <- if (is.finite(extent) && extent > 0) 0.04 * extent else 0.1
  threshold <- min(max(asNumericOr(threshold, CRAVE_COMMUNITY_MEMBERSHIP_DEFAULT), 0), 1)

  paths <- lapply(seq_len(ncol(membership)), function(j) {
    inside <- which(membership[, j] >= threshold)
    if (length(inside) == 0) return(NULL)
    path <- ellipsePath(layout$x[inside], layout$y[inside], level = level, pad = pad)
    if (is.null(path)) return(NULL)
    path$Community <- colnames(membership)[j]
    path$Members <- length(inside)
    path
  })
  paths <- Filter(Negate(is.null), paths)
  if (length(paths) == 0) return(NULL)
  out <- bind_rows(paths)
  out$Community <- factor(out$Community, levels = colnames(membership))
  out
}

#' The community network diagram.
#'
#' The same graph the Network tab draws when community detection is off -- nodes
#' and edges, laid out by Fruchterman-Reingold -- with each community ringed and
#' each gene coloured by the community it sits most in. Rendered through ggplot
#' and plotly rather than visNetwork, which has no way to draw a region around a
#' set of nodes; the trade is that the picture is not dragged about by the mouse,
#' and hovering still names the gene, its community and its membership.
#'
#' @param edges Both directions, as the adjacency list keeps them; deduplicated
#'   here so each pair is drawn, and laid out, once.
communityNetworkPlot <- function(edges, result, similarityName,
                                 threshold = CRAVE_COMMUNITY_MEMBERSHIP_DEFAULT,
                                 style = list()) {
  pairs <- edges %>%
    mutate(.a = pmin(from, to), .b = pmax(from, to)) %>%
    distinct(.a, .b, .keep_all = TRUE)

  nodes <- unique(c(pairs$.a, pairs$.b))
  g <- igraph::graph_from_data_frame(
    pairs %>% transmute(from = .a, to = .b, weight = abs(corr)),
    directed = FALSE, vertices = nodes)

  # Seeded through withSeed() so the same selection draws the same picture, and
  # the global stream is restored afterwards.
  coords <- withSeed(CRAVE_COMMUNITY_SEED, igraph::layout_with_fr(g))
  layout <- tibble(
    Gene = nodes,
    x = coords[, 1],
    y = coords[, 2],
    Community = factor(colnames(result$membership)[result$assignment[nodes]],
                       levels = colnames(result$membership)),
    Membership = round(unname(result$probability[nodes]), 3)
  )

  position <- setNames(seq_len(nrow(layout)), layout$Gene)
  segments <- tibble(
    x    = layout$x[position[pairs$.a]],
    y    = layout$y[position[pairs$.a]],
    xend = layout$x[position[pairs$.b]],
    yend = layout$y[position[pairs$.b]],
    corr = pairs$corr
  )

  ellipses <- communityEllipses(layout, result, threshold = threshold)

  p <- ggplot()
  if (!is.null(ellipses)) {
    # show.legend = FALSE, and the fill guide dropped below: the nodes and the
    # ellipses carry the same mapping, and a plot with both keys produces a
    # merged legend reading "(Community 1,1)" once it has been through plotly.
    p <- p + geom_polygon(
      data = ellipses,
      aes(x = x, y = y, group = Community, fill = Community, Members = Members),
      alpha = 0.15, colour = NA, show.legend = FALSE)
  }
  p <- p +
    geom_segment(data = segments,
                 aes(x = x, y = y, xend = xend, yend = yend, Similarity = corr),
                 colour = "#BBBBBB", linewidth = 0.3) +
    geom_point(data = layout,
               aes(x = x, y = y, colour = Community, Gene = Gene,
                   Membership = Membership),
               size = 3)

  # Labels only while they can be read; past that the hover does the work.
  if (nrow(layout) <= 60) {
    p <- p + geom_text(data = layout, aes(x = x, y = y, label = Gene),
                       vjust = -1, size = 3)
  }

  p <- p +
    theme_void() +
    guides(fill = "none") +
    # Two lines: plotly gives a title one line and clips it, and the second half
    # -- which measures produced this picture -- is the half a reader needs.
    labs(colour = "Community",
         title = communityTitle(result, similarityName, nrow(layout)))

  dims <- plotDims(style, CORRELATE_DEFAULTS)
  list(p = ggplotly(p, height = dims$height, width = dims$width),
       height = dims$height, layout = layout)
}

#' Attach community columns to the adjacency list.
#'
#' Both endpoints get their community, and under a soft method their membership as
#' well, so a row of the table says what it is joining and how firmly each end
#' belongs where it has been put.
withCommunityColumns <- function(out, result) {
  if (is.null(result) || !isTRUE(result$ok)) return(out)
  a <- as.character(out$`Gene A`)
  b <- as.character(out$`Gene B`)
  out$`Community (Gene A)` <- unname(result$assignment[a])
  out$`Community (Gene B)` <- unname(result$assignment[b])
  if (isTRUE(result$probabilistic)) {
    out$`Membership probability (Gene A)` <- round(unname(result$probability[a]), 4)
    out$`Membership probability (Gene B)` <- round(unname(result$probability[b]), 4)
  }
  out
}
