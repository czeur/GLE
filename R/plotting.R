# Plotting utilities

library(ggplot2)
library(igraph)

#' Set default graph visualization parameters
set_graph_parameters <- function(graph) {
  V(graph)$color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 0.4)
  V(graph)$frame.color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 1)
  V(graph)$label.color <- "black"
  V(graph)$size <- 15
  E(graph)$width <- 2
  E(graph)$color <- "darkgrey"
  graph
}

#' Save a ggplot with controlled panel size
save_myplot <- function(plt, plt_nm, width, height,
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)) dir.create(dir_name, recursive = TRUE)

  panel_plt <- egg::set_panel_size(
    p = plt, width = unit(width, "in"), height = unit(height, "in"))

  if (cairo) {
    ggsave(plt_nm, panel_plt,
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = "in",
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, panel_plt,
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = "in")
  }

  if (crop) knitr::plot_crop(plt_nm)
}

my_palette <- list(
  red = "#D55E00",
  blue = "#0072B2",
  green = "#009E73",
  yellow = "#E69F00",
  pink = "#CC79A7",
  light_blue = "#56B4E9",
  grey = "#999999",
  background = "#332288"
)

#' Set the default ggplot theme
theme_fct <- function(font_size1 = 11, font_size2 = 7.5) {
  theme_set(theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.caption = element_text(size = font_size2, hjust = 0,
                                  margin = margin(t = 15)),
      text = element_text(size = font_size1),
      axis.ticks = element_blank(),
      axis.text = element_text(size = font_size1),
      panel.grid.major = element_line(linewidth = 0.25)
    ))
}
