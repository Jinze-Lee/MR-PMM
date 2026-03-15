## MR-PMM TUTORIAL FUNCTIONS

#' Simulate multi-trait data under an MR-PMM style covariance structure.
#'
#' @param n Integer. Number of species / taxa.
#' @param m Integer. Number of response traits.
#' @param Sigma_phy Numeric matrix. Trait-level phylogenetic covariance matrix (m x m).
#' @param Sigma_res Numeric matrix. Trait-level residual covariance matrix (m x m).
#' @param phy A phylo object (ape) with n tips.
#'
#' @return A data.frame containing `animal` plus simulated traits `y1 ... ym`.
#'
#' @details
#' This helper function mirrors the generative logic of an MR-PMM:
#' 1) Build a phylogenetic covariance matrix C from the tree,
#' 2) Draw phylogenetic effects u ~ MVN(0, Sigma_phy ⊗ C),
#' 3) Draw independent residuals e ~ MVN(0, Sigma_res ⊗ I),
#' 4) Combine to produce each trait response.
mrpmm_sim <- function(n, m, Sigma_phy, Sigma_res, phy) {
  n <- n  # number of species
  m <- m  # number of traits

  # C: expected covariance among species due to shared ancestry.
  C <- ape::vcv.phylo(phy, corr = TRUE)

  # I: independent variance structure at the residual level.
  I <- diag(n)

  # Simulate phylogenetic random effects for all traits and all species.
  u <- MASS::mvrnorm(n = 1, mu = rep(0, n * m), Sigma = kronecker(Sigma_phy, C))

  # Simulate non-phylogenetic residual variation.
  e <- MASS::mvrnorm(n = 1, mu = rep(0, n * m), Sigma = kronecker(Sigma_res, I))

  # Number of species per trait block in the flattened vectors.
  rep_n <- length(u) / m

  # Create base output table indexed by species labels from the phylogeny.
  d <- data.frame(animal = phy$tip.label)

  # Split long vectors into trait-specific slices.
  u.list <- vector("list", m)
  e.list <- vector("list", m)

  for (i in seq_len(m)) {
    idx_start <- (1 + rep_n * i) - rep_n
    idx_end <- rep_n * i
    u.list[[i]] <- u[idx_start:idx_end]
    e.list[[i]] <- e[idx_start:idx_end]
  }

  names(u.list) <- paste0("u", seq_len(m))
  names(e.list) <- paste0("e", seq_len(m))

  # Assemble trait responses: y_i = u_i + e_i.
  for (i in seq_len(m)) {
    d[[paste0("y", i)]] <- unlist(u.list[i]) + unlist(e.list[i])
  }

  return(d)
}


#' Extract posterior correlation and partial-correlation matrices from MR-PMM VCV draws.
#'
#' @param vcv Posterior draws of covariance parameters (e.g., `fit$VCV`).
#' @param n_resp Integer. Number of response traits.
#' @param part.1 Character. Name of the phylogenetic random-effect block in `vcv`.
#' @param part.2 Character. Name of the residual/random-effect block in `vcv`.
#'
#' @return A list with two elements:
#' - `posteriors`: posterior draws for correlation and partial correlation.
#' - `matrices`: posterior-mean matrices at both levels.
#'
#' @details
#' This function is designed for interpretation workflows:
#' - `phy_cor` and `ind_cor` describe marginal trait associations,
#' - `par_phy_cor` and `par_ind_cor` describe conditional (partial) associations.
level_cor_vcv <- function(vcv, n_resp, part.1 = "animal", part.2 = "units") {

  # Index lower-triangle elements to keep unique trait-pair correlations.
  m <- matrix(seq_len(n_resp^2), n_resp, n_resp)
  for (i in seq_len(nrow(m))) {
    for (j in i:ncol(m)) {
      m[i, j] <- 0
    }
  }
  cols <- unique(as.vector(m))
  cols <- cols[cols != 0]

  # ---------- Phylogenetic level ----------
  draws <- vcv %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(tolower) %>%
    dplyr::rename_with(~ stringr::str_remove_all(., "trait"))

  cor_draws <- draws %>%
    dplyr::select(dplyr::contains(part.1)) %>%
    dplyr::rename_with(~ stringr::str_remove_all(., paste0(".", part.1)))

  phy_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  par_phy_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))

  for (i in seq_len(nrow(cor_draws))) {
    cor <- unlist(cor_draws[i, ])
    phy_res <- matrix(0, n_resp, n_resp)
    phy_res[] <- cor
    phy_res <- stats::cov2cor(phy_res)

    # Partial correlation: direct association between two traits
    # after conditioning on all other traits.
    par_phy_res <- corpcor::cor2pcor(phy_res)

    phy_cor[i, ] <- as.vector(phy_res)
    par_phy_cor[i, ] <- as.vector(par_phy_res)
  }

  phy_cor <- phy_cor %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(cols)) %>%
    setNames(names(cor_draws %>% dplyr::select(all_of(cols))))

  par_phy_cor <- par_phy_cor %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(cols)) %>%
    setNames(names(cor_draws %>% dplyr::select(all_of(cols))))

  phy_res[lower.tri(phy_res)] <- phy_cor %>%
    dplyr::summarise_all(~ mean(.)) %>%
    dplyr::slice(1) %>%
    as.numeric()

  phy_res[upper.tri(phy_res)] <- t(phy_res)[upper.tri(phy_res)]

  dimnames(phy_res)[[1]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))
  dimnames(phy_res)[[2]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))

  par_phy_res[lower.tri(par_phy_res)] <- par_phy_cor %>%
    dplyr::summarise_all(~ mean(.)) %>%
    dplyr::slice(1) %>%
    as.numeric()

  par_phy_res[upper.tri(par_phy_res)] <- t(par_phy_res)[upper.tri(par_phy_res)]

  dimnames(par_phy_res)[[1]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))
  dimnames(par_phy_res)[[2]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))

  # ---------- Residual / independent level ----------
  draws <- vcv %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(tolower) %>%
    dplyr::rename_with(~ stringr::str_remove_all(., "trait"))

  cor_draws <- draws %>%
    dplyr::select(dplyr::contains(part.2)) %>%
    dplyr::rename_with(~ stringr::str_remove_all(., paste0(".", part.2)))

  ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  par_ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))

  for (i in seq_len(nrow(cor_draws))) {
    cor <- unlist(cor_draws[i, ])
    ind_res <- matrix(0, n_resp, n_resp)
    ind_res[] <- cor
    ind_res <- stats::cov2cor(ind_res)

    par_ind_res <- corpcor::cor2pcor(ind_res)
    ind_cor[i, ] <- as.vector(ind_res)
    par_ind_cor[i, ] <- as.vector(par_ind_res)
  }

  ind_cor <- ind_cor %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(cols)) %>%
    setNames(names(cor_draws %>% dplyr::select(all_of(cols))))

  par_ind_cor <- par_ind_cor %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(cols)) %>%
    setNames(names(cor_draws %>% dplyr::select(all_of(cols))))

  ind_res[lower.tri(ind_res)] <- ind_cor %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.))) %>%
    dplyr::slice(1) %>%
    as.numeric()

  ind_res[upper.tri(ind_res)] <- t(ind_res)[upper.tri(ind_res)]

  dimnames(ind_res)[[1]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))
  dimnames(ind_res)[[2]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))

  par_ind_res[lower.tri(par_ind_res)] <- par_ind_cor %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.))) %>%
    dplyr::slice(1) %>%
    as.numeric()

  par_ind_res[upper.tri(par_ind_res)] <- t(par_ind_res)[upper.tri(par_ind_res)]

  dimnames(par_ind_res)[[1]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))
  dimnames(par_ind_res)[[2]] <- gsub(":", "", gsub("trait", "", gsub("[^:]+$", "", colnames(vcv)[1:n_resp])))

  result <- list(
    posteriors = list(
      phy_cor = phy_cor,
      par_phy_cor = par_phy_cor,
      ind_cor = ind_cor,
      par_ind_cor = par_ind_cor
    ),
    matrices = list(
      phy_mat = phy_res,
      par_phy_mat = par_phy_res,
      ind_mat = ind_res,
      par_ind_mat = par_ind_res
    )
  )

  return(result)
}
