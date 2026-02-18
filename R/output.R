# Experiment tracking: timestamped output, parameter logging

#' Create a timestamped output directory
#'
#' @param prefix Optional prefix for the directory name
#' @param base_dir Base output directory (default: "output")
#' @return Path to the created directory
create_output_dir <- function(prefix = NULL, base_dir = "output") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  dir_name <- if (!is.null(prefix)) {
    paste0(prefix, "_", timestamp)
  } else {
    timestamp
  }
  path <- file.path(base_dir, dir_name)
  dir.create(path, recursive = TRUE)
  path
}

#' Save experiment configuration as YAML
#'
#' @param config Named list of parameters
#' @param output_dir Directory to save to
save_config <- function(config, output_dir) {
  lines <- vapply(names(config), function(nm) {
    val <- config[[nm]]
    if (is.numeric(val) && length(val) > 1) {
      paste0(nm, ": [", paste(val, collapse = ", "), "]")
    } else {
      paste0(nm, ": ", val)
    }
  }, character(1))
  writeLines(lines, file.path(output_dir, "config.yaml"))
}

#' Save results in multiple formats
#'
#' @param results Object to save
#' @param name Base filename (without extension)
#' @param output_dir Directory to save to
#' @param formats Vector of formats: "rds", "rdata"
save_results <- function(results, name, output_dir, formats = "rds") {
  if ("rds" %in% formats) {
    saveRDS(results, file.path(output_dir, paste0(name, ".rds")))
  }
  if ("rdata" %in% formats) {
    save(results, file = file.path(output_dir, paste0(name, ".Rdata")))
  }
}

#' Save a plot in multiple formats
#'
#' @param plot_fn Function that produces the plot (called with no arguments)
#' @param name Base filename (without extension)
#' @param output_dir Directory to save to
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param formats Vector of formats: "pdf", "png"
save_plot <- function(plot_fn, name, output_dir,
                      width = 8, height = 5, formats = c("pdf", "png")) {
  if ("pdf" %in% formats) {
    pdf(file.path(output_dir, paste0(name, ".pdf")), width = width, height = height)
    plot_fn()
    dev.off()
  }
  if ("png" %in% formats) {
    png(file.path(output_dir, paste0(name, ".png")),
        width = width, height = height, units = "in", res = 150)
    plot_fn()
    dev.off()
  }
}
