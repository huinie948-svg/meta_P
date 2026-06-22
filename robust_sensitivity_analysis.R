library(metafor)
library(clubSandwich)

dir.create("robustness_results", showWarnings = FALSE)

optimizer_control <- list(
  optimizer = "optim",
  optmethod = "BFGS",
  maxit = 20000
)

datasets <- list(
  AP = "AP/ap_effects.csv",
  Phosphatase = "phosphatase/ap_effects.csv",
  TP = "TP/ap_effects.csv"
)
requested_endpoints <- Sys.getenv("ENDPOINTS", unset = "")
if (nzchar(requested_endpoints)) {
  keep <- trimws(strsplit(requested_endpoints, ",", fixed = TRUE)[[1]])
  datasets <- datasets[names(datasets) %in% keep]
}

read_endpoint <- function(path) {
  dat <- read.csv(path, fileEncoding = "latin1", check.names = TRUE)
  dat$yi_lnRR <- dat$RR
  dat$vi_lnRR <- dat$Vi
  dat$EffectID <- seq_len(nrow(dat))
  dat$study.id <- factor(dat$study.id)
  dat$EffectID <- factor(dat$EffectID)
  dat
}

extract_standard <- function(model, endpoint, analysis, studies) {
  pred <- predict(model)
  data.frame(
    endpoint = endpoint,
    analysis = analysis,
    k = model$k,
    studies = studies,
    estimate = unname(coef(model)[1]),
    se = model$se[1],
    df = NA_real_,
    p = model$pval[1],
    ci_lb = model$ci.lb[1],
    ci_ub = model$ci.ub[1],
    pi_lb = pred$pi.lb[1],
    pi_ub = pred$pi.ub[1],
    tau2_study = if (length(model$sigma2) >= 1) model$sigma2[1] else NA_real_,
    tau2_effect = if (length(model$sigma2) >= 2) model$sigma2[2] else NA_real_
  )
}

extract_robust <- function(model, robust_model, endpoint, analysis, studies) {
  pred <- predict(model)
  data.frame(
    endpoint = endpoint,
    analysis = analysis,
    k = model$k,
    studies = studies,
    estimate = unname(robust_model$beta[1]),
    se = robust_model$se[1],
    df = robust_model$dfs[1],
    p = robust_model$pval[1],
    ci_lb = robust_model$ci.lb[1],
    ci_ub = robust_model$ci.ub[1],
    pi_lb = pred$pi.lb[1],
    pi_ub = pred$pi.ub[1],
    tau2_study = model$sigma2[1],
    tau2_effect = model$sigma2[2]
  )
}

make_correlated_v <- function(vi, cluster, rho) {
  V <- matrix(0, nrow = length(vi), ncol = length(vi))
  groups <- split(seq_along(vi), cluster)
  for (idx in groups) {
    block <- rho * sqrt(outer(vi[idx], vi[idx]))
    diag(block) <- vi[idx]
    V[idx, idx] <- block
  }
  V
}

all_results <- list()
all_moderators <- list()

for (endpoint in names(datasets)) {
  message("Starting endpoint: ", endpoint)
  dat <- read_endpoint(datasets[[endpoint]])

  original <- rma.mv(
    yi = yi_lnRR,
    V = vi_lnRR,
    random = ~ 1 | study.id,
    data = dat,
    method = "REML",
    control = optimizer_control
  )
  all_results[[length(all_results) + 1]] <- extract_standard(
    original, endpoint, "Original study-random-intercept model",
    length(unique(dat$study.id))
  )

  multilevel <- rma.mv(
    yi = yi_lnRR,
    V = vi_lnRR,
    random = ~ 1 | study.id/EffectID,
    data = dat,
    method = "REML",
    control = optimizer_control
  )
  all_results[[length(all_results) + 1]] <- extract_standard(
    multilevel, endpoint, "Two-level study/effect random model",
    length(unique(dat$study.id))
  )

  robust_fit <- robust(
    multilevel,
    cluster = dat$study.id,
    clubSandwich = TRUE
  )
  all_results[[length(all_results) + 1]] <- extract_robust(
    multilevel, robust_fit, endpoint,
    "Two-level model with StudyID-clustered CR2 inference",
    length(unique(dat$study.id))
  )

  complete_var <- dat[!is.na(dat$sc) & !is.na(dat$se), ]
  complete_model <- rma.mv(
    yi = yi_lnRR,
    V = vi_lnRR,
    random = ~ 1 | study.id/EffectID,
    data = complete_var,
    method = "REML",
    control = optimizer_control
  )
  complete_robust <- robust(
    complete_model,
    cluster = complete_var$study.id,
    clubSandwich = TRUE
  )
  all_results[[length(all_results) + 1]] <- extract_robust(
    complete_model,
    complete_robust,
    endpoint,
    "CR2 model excluding observations with missing SD/SE",
    length(unique(complete_var$study.id))
  )

  for (rho in c(0.3, 0.6)) {
    message("  Correlated sampling model rho=", rho)
    tryCatch({
      V_rho <- make_correlated_v(dat$vi_lnRR, dat$study.id, rho)
      rho_model <- rma.mv(
        yi = yi_lnRR,
        V = V_rho,
        random = ~ 1 | study.id/EffectID,
        data = dat,
        method = "REML",
        sparse = TRUE,
        control = optimizer_control
      )
      rho_robust <- robust(
        rho_model,
        cluster = dat$study.id,
        clubSandwich = TRUE
      )
      all_results[[length(all_results) + 1]] <- extract_robust(
        rho_model,
        rho_robust,
        endpoint,
        paste0("CR2 model with assumed within-study sampling correlation rho=", rho),
        length(unique(dat$study.id))
      )
    }, error = function(e) {
      warning(
        sprintf(
          "Skipping %s rho=%s model: %s",
          endpoint, rho, conditionMessage(e)
        )
      )
    })
  }

  moderator_names <- c(
    "inoculant.type", "experimental.type", "stress",
    "fertilizer", "extraction.part"
  )
  for (moderator in moderator_names) {
    use <- dat[!is.na(dat[[moderator]]) & dat[[moderator]] != "unknow", ]
    use[[moderator]] <- droplevels(factor(use[[moderator]]))
    if (nlevels(use[[moderator]]) < 2) next
    form <- as.formula(paste0("~ 0 + ", moderator))
    tryCatch({
      mod_fit <- rma.mv(
        yi = yi_lnRR,
        V = vi_lnRR,
        mods = form,
        random = ~ 1 | study.id/EffectID,
        data = use,
        method = "REML",
        control = optimizer_control
      )
      mod_robust <- robust(
        mod_fit,
        cluster = use$study.id,
        clubSandwich = TRUE
      )
      coefficient_names <- rownames(mod_robust)
      if (is.null(coefficient_names)) {
        coefficient_names <- names(coef(mod_fit))
      }
      all_moderators[[length(all_moderators) + 1]] <- data.frame(
        endpoint = endpoint,
        moderator = moderator,
        level = coefficient_names,
        estimate = unname(mod_robust$beta),
        se = mod_robust$se,
        df = mod_robust$dfs,
        p = mod_robust$pval,
        ci_lb = mod_robust$ci.lb,
        ci_ub = mod_robust$ci.ub,
        k = nrow(use),
        studies = length(unique(use$study.id))
      )
    }, error = function(e) {
      warning(
        sprintf(
          "Skipping %s moderator %s: %s",
          endpoint, moderator, conditionMessage(e)
        )
      )
    })
  }

  write.csv(
    do.call(rbind, all_results),
    "robustness_results/overall_sensitivity_results_checkpoint.csv",
    row.names = FALSE
  )
  if (length(all_moderators) > 0) {
    write.csv(
      do.call(rbind, all_moderators),
      "robustness_results/moderator_CR2_results_checkpoint.csv",
      row.names = FALSE
    )
  }
  message("Completed endpoint: ", endpoint)
}

results <- do.call(rbind, all_results)
moderators <- do.call(rbind, all_moderators)

write.csv(
  results,
  "robustness_results/overall_sensitivity_results.csv",
  row.names = FALSE
)
write.csv(
  moderators,
  "robustness_results/moderator_CR2_results.csv",
  row.names = FALSE
)

sink("robustness_results/session_info.txt")
print(sessionInfo())
sink()

print(results)
