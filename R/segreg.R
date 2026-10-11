#' Segmented regression with automatic selection of the number of breakpoints
#'
#' Fits models with 0 to `k_max` breakpoints, selects the best one by AIC or
#' BIC, and assigns each observation to a phase.
#'
#' @param formula A formula of the form `y ~ x` (one continuous predictor).
#' @param data A data.frame.
#' @param k_max Maximum number of breakpoints to test. The linear model
#'   (0 breakpoints) is always included in the comparison. If the selected
#'   number of breakpoints equals `k_max`, a warning is issued because more
#'   breakpoints may exist; increase `k_max` in that case.
#' @param criterion Model selection criterion: "AIC" or "BIC".
#' @param n_grid Number of grid points used to draw the fitted line.
#' @param ... Additional arguments passed to `segmented::segmented()`
#'   (e.g., `control = seg.control(...)`).
#' @return An object of class "segreg".
#' @examples
#' \donttest{
#' # data upload
#' if(!require(readr)) install.packages("readr")
#' library (readr)
#' df= read.csv (paste0("https://raw.githubusercontent.com/agronomy4future/",
#'     "raw_data_practice/refs/heads/main/Practice_Segmented_Regression.csv"))
#' Mean= aggregate (yield ~ year, data= df, FUN= mean)
#' names(Mean)[2]= "yield_Mean"
#'
#' # import library
#' if(!require(remotes)) install.packages("remotes")
#' if (!requireNamespace("segreg", quietly = TRUE)) {
#'  remotes::install_github("agronomy4future/segreg", force= TRUE)
#' }
#' library(remotes)
#' library(segreg)
#'
#' # run the model
#' model= segreg (yield_Mean ~ year, data= Mean)
#' print (model) # Selected breakpoints, model comparison table, and slope by phase
#' summary (model) # Full model output with breakpoint CIs and slope CIs
#' seg_phase (model) # Original data with an added phase column
#'
#' if(!require(ggplot2)) install.packages("ggplot2")
#' library(ggplot2)
#'
#' # ggplot
#' ggplot(data= seg_phase(model), aes (x= year, y= yield_Mean)) +
#'   geom_bar(aes(fill=phase), stat="identity", position="dodge", width=0.9) +
#'   segreg_line(model, linetype= "dashed", linewidth= 1) +
#'   segreg_vline(model, type= "phase", linewidth= 1)
#' }
#' # Github: https://github.com/agronomy4future/segreg
#' # All Rights Reserved (c) J.K Kim (kimjk@agronomy4future.com)
#' @export
segreg = function(formula, data, k_max = 3, criterion = c("AIC", "BIC"),
                  n_grid = 200, ...) {
  data_expr = substitute(data)          # for display: data = Mean
  criterion = match.arg(criterion)
  vars = all.vars(formula)
  if (length(vars) != 2) stop("formula must be of the form y ~ x")
  yname = vars[1]; xname = vars[2]

  data = as.data.frame(data)
  data = data[stats::complete.cases(data[, c(yname, xname)]), , drop = FALSE]

  # Embed the data in the call so segmented() can find it inside this function
  lin = do.call("lm", list(formula = formula, data = data))
  seg_z = stats::as.formula(paste("~", xname))
  ic = if (criterion == "AIC") stats::AIC else stats::BIC

  fit_k = function(k, ...) {
    s = tryCatch(
      suppressWarnings(segmented::segmented(lin, seg.Z = seg_z, npsi = k, ...)),
      error = function(e) NULL)
    if (inherits(s, "segmented") && !is.null(s$psi) && nrow(s$psi) == k) s else NULL
  }
  dots = list(...)
  fits = list(`0` = lin)
  for (k in seq_len(k_max)) {
    s = fit_k(k, ...)
    # Bootstrap restarts can fail with small samples; retry with n.boot = 0
    if (is.null(s) && is.null(dots$control))
      s = fit_k(k, control = segmented::seg.control(n.boot = 0))
    if (!is.null(s)) fits[[as.character(k)]] = s
  }

  ks = 0:k_max
  ic_table = data.frame(k = ks, value = vapply(ks, function(k) {
    f = fits[[as.character(k)]]
    if (is.null(f)) NA_real_ else ic(f)
  }, numeric(1)))
  names(ic_table)[2] = criterion

  best_k = ic_table$k[which.min(ic_table[[criterion]])]
  best = fits[[as.character(best_k)]]
  bp = if (best_k > 0) sort(unname(best$psi[, "Est."])) else numeric(0)

  # Warn when the selection hits the upper limit: more breakpoints may exist
  at_kmax = k_max > 0 && best_k == k_max
  if (at_kmax)
    warning(sprintf(paste0("The selected number of breakpoints equals k_max (%d). ",
                           "More breakpoints may exist; consider increasing k_max ",
                           "(e.g., k_max = %d)."), k_max, k_max + 2), call. = FALSE)

  data$phase = cut(data[[xname]], breaks = c(-Inf, bp, Inf),
                   labels = paste0("P", seq_len(length(bp) + 1)))

  # Include breakpoints in the grid so the kinks are drawn exactly
  x = data[[xname]]
  xg = sort(unique(c(seq(min(x), max(x), length.out = n_grid), bp)))
  line = data.frame(xg)
  names(line) = xname
  line$fit = as.numeric(stats::predict(best, newdata = line))

  # Clean up calls for display: lm(formula = y ~ x, data = Mean)
  lm_call = bquote(lm(formula = .(formula), data = .(data_expr)))
  fits[["0"]]$call = lm_call
  for (k in setdiff(names(fits), "0"))
    fits[[k]]$call = bquote(segmented(.(lm_call), seg.Z = .(seg_z), npsi = .(as.numeric(k))))
  best = fits[[as.character(best_k)]]

  structure(list(model = best, models = fits, best_k = best_k,
                 criterion = criterion, ic_table = ic_table,
                 breakpoints = bp, at_kmax = at_kmax, data = data, line = line,
                 xname = xname, yname = yname, call = match.call()),
            class = "segreg")
}

#' @export
print.segreg = function(x, ...) {
  fmt = function(v) vapply(round(v, 2), function(z) format(z), "")
  bp = x$breakpoints
  cat("Segmented regression: ", x$yname, " ~ ", x$xname, "   (",
      x$criterion, ", ", x$best_k, " breakpoint", if (x$best_k != 1) "s", ")\n\n",
      sep = "")

  # Model comparison table (mark the selected k)
  tab = x$ic_table
  tab[[x$criterion]] = ifelse(is.na(tab[[x$criterion]]), "fail",
                              formatC(tab[[x$criterion]], format = "f", digits = 4))
  tab$` ` = ifelse(tab$k == x$best_k, "<- selected", "")
  print(tab, row.names = FALSE, right = TRUE)

  cat("\nBreakpoints:", if (length(bp)) paste(fmt(bp), collapse = ", ") else "none", "\n")
  if (isTRUE(x$at_kmax))
    cat("Note: the selected k equals k_max (", x$best_k,
        "); more breakpoints may exist. Try a larger k_max.\n", sep = "")

  # Slope by phase
  xr = range(x$data[[x$xname]])
  edges = c(xr[1], bp, xr[2])
  slopes = if (x$best_k > 0) unname(segmented::slope(x$model)[[1]][, "Est."])
  else unname(stats::coef(x$model)[x$xname])
  sl = data.frame(phase = paste0("P", seq_along(slopes)),
                  range = paste(fmt(utils::head(edges, -1)), fmt(edges[-1]), sep = "-"),
                  slope = round(slopes, 4))
  cat("\nSlope by phase:\n")
  print(sl, row.names = FALSE)

  cat("\n")
  print(x$data[, c(x$xname, x$yname, "phase")])
  invisible(x)
}

#' @export
summary.segreg = function(object, ...) {
  cat("Selected breakpoints (", object$criterion, "): ", object$best_k, "\n", sep = "")
  print(object$ic_table, row.names = FALSE)
  cat("\n"); print(summary(object$model))
  if (object$best_k > 0) {
    cat("\nBreakpoint CI:\n"); print(stats::confint(object$model))
    cat("\nSlope by phase:\n"); print(segmented::slope(object$model))
  }
  invisible(object)
}

#' Fitted line data (x, fit) for custom plotting or export
#' @param model A segreg object.
#' @export
segreg_fit = function(model) {
  stopifnot(inherits(model, "segreg"))
  model$line
}

#' Estimated breakpoints
#' @param model A segreg object.
#' @export
seg_breaks = function(model) model$breakpoints

#' Original data with a phase column
#' @param model A segreg object.
#' @export
seg_phase = function(model) model$data

#' @exportS3Method ggplot2::fortify
fortify.segreg = function(model, data, ...) model$line

# Normalize color/col to colour and merge with defaults
.seg_args = function(args, defaults) {
  nm = names(args)
  if (!is.null(nm)) names(args)[nm %in% c("color", "col")] = "colour"
  utils::modifyList(defaults, args)
}

#' Fitted segmented line layer
#' @param model A segreg object.
#' @param ... Passed to `geom_line()` (color, linetype, linewidth, alpha, etc.).
#' @importFrom ggplot2 .data
#' @export
segreg_line = function(model, ...) {
  stopifnot(inherits(model, "segreg"))
  d = data.frame(x = model$line[[model$xname]], y = model$line$fit)
  args = .seg_args(list(...), list(colour = "black", linetype = "dashed", linewidth = 1))
  do.call(ggplot2::geom_line,
          c(list(data = d, mapping = ggplot2::aes(x = .data$x, y = .data$y),
                 inherit.aes = FALSE), args))
}

#' Phase boundaries: midpoints between the last and first observations of
#' adjacent phases (useful for separating groups in bar charts)
#' @param model A segreg object.
#' @export
seg_bounds = function(model) {
  stopifnot(inherits(model, "segreg"))
  d = model$data[order(model$data[[model$xname]]), ]
  x = d[[model$xname]]; ph = as.character(d$phase)
  idx = which(ph[-1] != ph[-length(ph)])
  (x[idx] + x[idx + 1]) / 2
}

#' Vertical line layer
#' @param model A segreg object.
#' @param type "breakpoint" = estimated breakpoint positions,
#'   "phase" = boundaries between phase groups.
#' @param ... Passed to `geom_vline()` (color, linetype, linewidth, alpha, etc.).
#' @export
segreg_vline = function(model, type = c("breakpoint", "phase"), ...) {
  stopifnot(inherits(model, "segreg"))
  type = match.arg(type)
  xi = if (type == "phase") seg_bounds(model) else model$breakpoints
  if (length(xi) == 0) return(NULL)
  args = .seg_args(list(...), list(colour = "steelblue", linetype = "dashed"))
  do.call(ggplot2::geom_vline, c(list(xintercept = xi), args))
}
