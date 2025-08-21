
fit_cohort = function(d,
                      fit_breakpoints = T,
                      with_initiation = T,
                      max_segments = 3,
                      min_segment_size = 3,
                      comparison = "bic",
                      enforce_rho_separation = T,
                      alpha_rho = .05,
                      models_to_fit = c("exponential", "logistic", "gompertz", "monomolecular", "quadraticexp"),
                      chains = 4,
                      iter = 1000,
                      seed = 123,
                      cores = 4,
                      t_prior_sd = 1,
                      method = "sampling",
                      use_elbo = F) {

  all_data = readRDS("/Users/jovoni/Dropbox/Zenodo_biPOD_v2/Zenodo/method_validation/sauer/data/processed_data.rds")
  data = all_data %>%
    dplyr::filter(OV04_patient == 828) %>%
    dplyr::filter(time <= 20, time >= -10)

  data$ID = data$mouse

  d = data
  unique_ids = unique(data$ID)
  id = unique_ids[1]

  cohort_fits = lapply(unique_ids, function(id) {
    d_id = d %>% dplyr::filter(ID == id) %>%
      dplyr::filter(count > 0)

    if (fit_breakpoints) {
      x_bp = biPOD::fit_breakpoints(
        data = d_id,
        with_initiation = with_initiation,
        max_segments = max_segments,
        min_segment_size = min_segment_size,
        comparison = comparison,
        enforce_rho_separation = enforce_rho_separation,
        alpha_rho = alpha_rho,
        models_to_fit = c("exponential"),
        chains = chains,
        iter = iter,
        seed = seed,
        cores = cores,
        t_prior_sd = t_prior_sd
      )

      breakpoints = x_bp$final_breakpoints
    } else {
      breakpoints = NULL
      x_bp = NULL
    }

    x = biPOD::fit_growth(
      data = d_id, breakpoints = breakpoints,
      with_initiation = with_initiation, chains = chains,
      iter = iter, seed = seed, cores = cores, comparison = comparison,
      models_to_fit = models_to_fit, method = method, use_elbo = use_elbo
    )

    fit = list(breakpoints_fit = x_bp, growth_fit = x, data = d_id)
    fit
  })

  names(cohort_fits) = unique_ids

  # Summarize cohort
  param_df_cohort = dplyr::tibble()
  breakpoints_df_cohort = dplyr::tibble()


  names(cohort_fits)

  id = 529
  for (id in unique_ids) {
    id = as.character(id)
    f = cohort_fits[[id]]

    param_df = f$growth_fit$fit$summary %>%
      dplyr::filter(!grepl("log_lik", variable)) %>%
      dplyr::filter(!grepl("yrep", variable)) %>%
      dplyr::filter(!grepl("mu_pred", variable)) %>%
      dplyr::filter(variable != "lp__") %>%
      dplyr::select(variable, median, sd) %>%
      dplyr::mutate(ID = id) %>%
      dplyr::mutate(growth = f$growth_fit$best_model)

    param_df_cohort = dplyr::bind_rows(param_df_cohort, param_df)

    if (fit_breakpoints) {
      breakpoints = sort(f$growth_fit$breakpoints)
      if (length(breakpoints) > 0) {
        breakpoints_df = dplyr::tibble(ID = id, idx = 1:length(breakpoints), breakpoints = breakpoints)
        breakpoints_df_cohort = dplyr::bind_rows(breakpoints_df_cohort, breakpoints_df)
      }
    }
  }

  x_min = min(d$time)
  x_max = max(d$time)

  d %>%
    ggplot(mapping = aes(x = time, y = time)) +
    geom_density2d()

  length = 500
  x_min = min(d$time)
  x_max = max(d$time)
  time_grid = seq(x_min, x_max, length = length)

  patient_curves = lapply(unique_ids, function(id) {
    id = as.character(id)
    f = cohort_fits[[id]]

    biPOD:::get_data_for_growth_plot(f$growth_fit, f$data, time_grid = time_grid) %>%
      dplyr::mutate(ID = id)
  }) %>% do.call("bind_rows", .)

}

cluster_cohort = function() {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(factoextra)
  library(dtw)  # for dynamic time warping
  library(TSclust)  # for time series clustering


  # ============================================================================
  # METHOD 1: FEATURE-BASED CLUSTERING
  # Extract summary features from each patient's trajectory
  # ============================================================================

  extract_trajectory_features <- function(data) {
    data %>%
      group_by(ID) %>%
      summarise(
        # Central tendency features
        mean_median = mean(median, na.rm = TRUE),
        median_median = median(median, na.rm = TRUE),
        final_value = last(median),
        initial_value = first(median),

        # Variability features
        sd_median = sd(median, na.rm = TRUE),
        cv_median = sd_median / abs(mean_median),  # coefficient of variation
        range_median = max(median, na.rm = TRUE) - min(median, na.rm = TRUE),

        # Trend features
        linear_slope = {
          if(n() > 2) {
            lm_fit <- lm(median ~ time)
            coef(lm_fit)[2]  # slope
          } else NA
        },

        # Uncertainty features
        mean_uncertainty = mean(upper - lower, na.rm = TRUE),
        uncertainty_trend = {
          if(n() > 2) {
            uncertainty <- upper - lower
            lm_fit <- lm(uncertainty ~ time)
            coef(lm_fit)[2]  # slope of uncertainty
          } else NA
        },

        # Shape features
        peak_time = time[median == max(median)],

        # Area under curve
        auc = sum(median * c(diff(time), 0), na.rm = TRUE),

        .groups = 'drop'
      ) %>%
      # Handle any remaining NAs
      mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
  }

  print("=== METHOD 1: FEATURE-BASED CLUSTERING ===")

  # Extract features
  dynamics_data = patient_curves
  feature_df <- extract_trajectory_features(dynamics_data)
  feature_df = feature_df %>% dplyr::select(ID, auc)

  print("Extracted features:")
  print(colnames(feature_df))

  # Prepare features for clustering (exclude ID)
  cluster_features <- feature_df %>%
    dplyr::select(-ID) %>%
    scale()

  # Determine optimal number of clusters
  print("Determining optimal number of clusters...")

  # Method A: Elbow method
  wss <- sapply(1:(nrow(cluster_features)-1), function(k) {
    kmeans(cluster_features, centers = k, nstart = 100)$tot.withinss
  })

  # Method B: Silhouette method
  sil_scores <- sapply(2:(nrow(cluster_features)-1), function(k) {
    km_fit <- kmeans(cluster_features, centers = k, nstart = 50)
    sil <- silhouette(km_fit$cluster, dist(cluster_features))
    mean(sil[,3])
  })

  optimal_k <- which.max(sil_scores) + 1
  print(paste("Optimal number of clusters (silhouette method):", optimal_k))

  # Perform K-means clustering
  set.seed(123)
  kmeans_fit <- kmeans(cluster_features, centers = optimal_k, nstart = 50)

  # Add cluster assignments
  feature_df$cluster <- kmeans_fit$cluster

  feature_df %>% dplyr::select(ID, cluster) %>%
    dplyr::left_join(data, by = "ID") %>%
    ggplot(mapping = aes(x = time, y = count, col = ID)) +
    geom_line() +
    facet_grid(group~cluster) +
    theme_bw() +
    geom_vline(xintercept = 10)

  plot_growth_fit(cohort_fits[[id]]$growth_fit, cohort_fits[[id]]$data)

  print("Feature-based clustering results:")
  print(table(feature_df$cluster))

  # ============================================================================
  # METHOD 2: TRAJECTORY SHAPE CLUSTERING (Dynamic Time Warping)
  # ============================================================================

  print("\n=== METHOD 2: TRAJECTORY SHAPE CLUSTERING ===")

  # Prepare data in wide format for DTW
  trajectory_matrix <- dynamics_data %>%
    select(ID, time, median) %>%
    pivot_wider(names_from = time, values_from = median, names_prefix = "t_") %>%
    column_to_rownames("ID") %>%
    as.matrix()

  # Handle missing values
  trajectory_matrix[is.na(trajectory_matrix)] <- 0

  # Calculate DTW distance matrix
  print("Calculating DTW distances...")

  # Method A: Using dtw package directly
  library(dtw)
  n_patients <- nrow(trajectory_matrix)
  patient_names <- rownames(trajectory_matrix)

  # Create DTW distance matrix
  dtw_dist_matrix <- matrix(0, n_patients, n_patients)
  rownames(dtw_dist_matrix) <- patient_names
  colnames(dtw_dist_matrix) <- patient_names

  for(i in 1:(n_patients-1)) {
    for(j in (i+1):n_patients) {
      dtw_result <- dtw(trajectory_matrix[i,], trajectory_matrix[j,],
                        distance.only = TRUE)
      dtw_dist_matrix[i,j] <- dtw_result$distance
      dtw_dist_matrix[j,i] <- dtw_result$distance
    }
    if(i %% 10 == 0) print(paste("Processed", i, "of", n_patients-1, "patients"))
  }

  # Convert to dist object
  dtw_distances <- as.dist(dtw_dist_matrix)

  # Alternative Method B: Using TSclust package (faster for large datasets)
  # library(TSclust)
  # dtw_distances <- diss(trajectory_matrix, "DTWARP")

  # Alternative faster DTW clustering using TSclust
  print("\n=== ALTERNATIVE: FASTER DTW WITH TSclust ===")

  tryCatch({
    library(TSclust)

    # TSclust method - much faster for large datasets
    print("Using TSclust for faster DTW computation...")
    dtw_distances_fast <- diss(trajectory_matrix, "DTWARP")

    # Hierarchical clustering
    hclust_dtw_fast <- hclust(dtw_distances_fast, method = "ward.D2")
    dtw_clusters_fast <- cutree(hclust_dtw_fast, k = optimal_k)

    print("Fast DTW clustering results:")
    print(table(dtw_clusters_fast))

    # Use the faster method if available
    dtw_clusters <- dtw_clusters_fast
    hclust_dtw <- hclust_dtw_fast

  }, error = function(e) {
    print("TSclust package not available, using manual DTW calculation")
    print("Install with: install.packages('TSclust')")
  })

  # ============================================================================
  # METHOD 3: FUNCTIONAL DATA CLUSTERING
  # ============================================================================

  print("\n=== METHOD 3: FUNCTIONAL DATA CLUSTERING ===")

  # Smooth trajectories and cluster based on functional form
  tryCatch({
    library(fda)
    library(funHDDC)

    # Create basis functions (B-splines)
    time_range <- range(dynamics_data$time)
    basis <- create.bspline.basis(time_range, nbasis = 7)

    # Convert to functional data
    trajectory_list <- dynamics_data %>%
      group_by(ID) %>%
      arrange(time) %>%
      group_split()

    # Smooth each trajectory
    fd_objects <- lapply(trajectory_list, function(patient_data) {
      smooth.basis(patient_data$time, patient_data$median, basis)$fd
    })

    # Combine into functional data object
    fd_combined <- do.call(cbind, lapply(fd_objects, function(x) x$coefs))
    fd_final <- fd(fd_combined, basis)

    # Functional clustering
    fclust_result <- funHDDC(fd_final, K = 2:optimal_k, model = "AkjBkQkDk")

    print("Functional clustering results:")
    print(table(fclust_result$class))

  }, error = function(e) {
    print("Functional clustering requires 'fda' and 'funHDDC' packages")
    print("Install with: install.packages(c('fda', 'funHDDC'))")
  })

  # ============================================================================
  # METHOD 4: MIXED EFFECTS CLUSTERING
  # ============================================================================

  print("\n=== METHOD 4: TRAJECTORY PARAMETERS CLUSTERING ===")

  # Fit individual models and cluster based on parameters
  individual_params <- dynamics_data %>%
    group_by(ID) %>%
    do(model = {
      if(nrow(.) > 3) {
        # Fit quadratic model to capture curvature
        lm(median ~ poly(time, 2), data = .)
      } else {
        # Fallback to linear if insufficient data
        lm(median ~ time, data = .)
      }
    }) %>%
    mutate(
      intercept = map_dbl(model, ~coef(.)[1]),
      linear_term = map_dbl(model, ~ifelse(length(coef(.)) > 1, coef(.)[2], 0)),
      quadratic_term = map_dbl(model, ~ifelse(length(coef(.)) > 2, coef(.)[3], 0)),
      r_squared = map_dbl(model, ~summary(.)$r.squared)
    ) %>%
    select(-model)

  # Cluster based on trajectory parameters
  param_features <- individual_params %>%
    select(-ID) %>%
    scale()

  param_clusters <- kmeans(param_features, centers = optimal_k, nstart = 50)
  individual_params$cluster <- param_clusters$cluster

  print("Parameter-based clustering results:")
  print(table(individual_params$cluster))

  # ============================================================================
  # VISUALIZATION AND COMPARISON
  # ============================================================================

  print("\n=== CREATING VISUALIZATIONS ===")

  # Combine all clustering results
  all_clusters <- tibble(ID = names(dtw_clusters), dtw_cluster = dtw_clusters)

  # Add clusters back to original data for plotting
  plot_data <- dynamics_data %>%
    left_join(all_clusters, by = "ID")

  # Function to create trajectory plots by cluster
  plot_trajectories <- function(data, cluster_col, title) {
    ggplot(data, aes(x = time, y = median, group = ID)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(!!sym(cluster_col))),
                  alpha = 0.1) +
      geom_line(aes(color = factor(!!sym(cluster_col))), alpha = 0.6) +
      facet_wrap(~factor(!!sym(cluster_col)), labeller = label_both) +
      labs(title = title,
           x = "Time",
           y = "Median Value",
           color = "Cluster",
           fill = "Cluster") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  # Create plots for each method
  if(exists("plot_data") && nrow(plot_data) > 0) {

    # Feature-based clustering plot
    p1 <- plot_trajectories(plot_data, "feature_cluster", "Feature-Based Clustering")
    print(p1)

    # DTW clustering plot
    p2 <- plot_trajectories(plot_data, "dtw_cluster", "DTW-Based Clustering")
    print(p2)

    # Parameter clustering plot
    p3 <- plot_trajectories(plot_data, "param_cluster", "Parameter-Based Clustering")
    print(p3)

    # Cluster comparison
    print("\n=== CLUSTER COMPARISON ===")
    cluster_agreement <- all_clusters %>%
      select(-ID) %>%
      cor(use = "complete.obs")
    print("Correlation between clustering methods:")
    print(round(cluster_agreement, 3))
  }

  # ============================================================================
  # RECOMMENDATIONS BASED ON YOUR DATA CHARACTERISTICS
  # ============================================================================

  print("\n=== CLUSTERING METHOD RECOMMENDATIONS ===")
  print("Choose based on your research question:")
  print("")
  print("1. FEATURE-BASED: Use when you care about overall trajectory characteristics")
  print("   - Good for: identifying 'high vs low responders', 'fast vs slow changers'")
  print("   - Pros: Interpretable features, handles noise well")
  print("   - Cons: May miss complex temporal patterns")
  print("")
  print("2. DTW-BASED: Use when trajectory SHAPE matters most")
  print("   - Good for: finding patients with similar temporal patterns")
  print("   - Pros: Captures shape similarity, handles time shifts")
  print("   - Cons: Computationally expensive, sensitive to noise")
  print("")
  print("3. PARAMETER-BASED: Use when you want to understand trajectory types")
  print("   - Good for: linear vs quadratic vs exponential patterns")
  print("   - Pros: Model-based interpretation, handles irregular sampling")
  print("   - Cons: Assumes specific functional forms")
  print("")
  print("4. FUNCTIONAL: Use for smooth, densely sampled trajectories")
  print("   - Good for: complex smooth patterns, derivative-based features")
  print("   - Pros: Mathematically rigorous, captures smooth patterns")
  print("   - Cons: Requires smooth data, complex interpretation")

  # Performance evaluation function
  evaluate_clustering <- function(clusters, data) {
    # Silhouette analysis
    sil <- silhouette(clusters, dist(data))
    avg_sil <- mean(sil[,3])

    # Within-cluster sum of squares
    wss <- sum(sapply(unique(clusters), function(k) {
      cluster_data <- data[clusters == k, , drop = FALSE]
      if(nrow(cluster_data) > 1) {
        sum(dist(cluster_data)^2) / (2 * nrow(cluster_data))
      } else 0
    }))

    list(silhouette = avg_sil, wss = wss)
  }

  print("\n=== CLUSTERING QUALITY METRICS ===")
  if(exists("feature_df") && exists("cluster_features")) {
    metrics1 <- evaluate_clustering(feature_df$cluster, cluster_features)
    print(paste("Feature-based - Silhouette:", round(metrics1$silhouette, 3)))

    metrics2 <- evaluate_clustering(dtw_clusters, trajectory_matrix)
    print(paste("DTW-based - Silhouette:", round(metrics2$silhouette, 3)))

    metrics3 <- evaluate_clustering(individual_params$cluster, param_features)
    print(paste("Parameter-based - Silhouette:", round(metrics3$silhouette, 3)))
  }
}




