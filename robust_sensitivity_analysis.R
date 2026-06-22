library(metafor)

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
  dat$study.id <- factor(dat$study.id)
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
    tau2_study = if (length(model$sigma2) >= 1) model$sigma2[1] else NA_real_
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

  complete_var <- dat[!is.na(dat$sc) & !is.na(dat$se), ]
  complete_model <- rma.mv(
    yi = yi_lnRR,
    V = vi_lnRR,
    random = ~ 1 | study.id,
    data = complete_var,
    method = "REML",
    control = optimizer_control
  )
  all_results[[length(all_results) + 1]] <- extract_standard(
    complete_model, endpoint,
    "StudyID model excluding observations with missing SD/SE",
    length(unique(complete_var$study.id))
  )

  for (rho in c(0.3, 0.6)) {
    message("  Correlated sampling model rho=", rho)
    tryCatch({
      V_rho <- make_correlated_v(dat$vi_lnRR, dat$study.id, rho)
      rho_model <- rma.mv(
        yi = yi_lnRR,
        V = V_rho,
        random = ~ 1 | study.id,
        data = dat,
        method = "REML",
        sparse = TRUE,
        control = optimizer_control
      )
      all_results[[length(all_results) + 1]] <- extract_standard(
        rho_model, endpoint,
        paste0("StudyID model with assumed within-study sampling correlation rho=", rho),
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

  write.csv(
    do.call(rbind, all_results),
    "robustness_results/overall_sensitivity_results_checkpoint.csv",
    row.names = FALSE
  )
  message("Completed endpoint: ", endpoint)
}

results <- do.call(rbind, all_results)

write.csv(
  results,
  "robustness_results/overall_sensitivity_results.csv",
  row.names = FALSE
)

sink("robustness_results/session_info.txt")
print(sessionInfo())
sink()

print(results)
