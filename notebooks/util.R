
######### Custom ggplot2 theme #########

font_size <- 8
font_size_small <- 6
font_size_tiny <- 5
font_size_large <- 10
publ_theme <- cowplot::theme_cowplot(font_size = font_size, rel_small = font_size_small / font_size,
                                     rel_tiny = font_size_tiny / font_size, rel_large = font_size_large / font_size,
                                     line_size = 0.3) +
  theme(plot.title = element_text(size = font_size),
        axis.title = element_text(size = font_size_small),
        legend.title = element_text(size = font_size_small),
        strip.background = element_blank(),
        strip.text = element_text(size = font_size_small), 
        panel.spacing.x = unit(2, "mm"))
theme_set(publ_theme)

tikzDevice::setTikzDefaults()
options(tikzDocumentDeclaration = union(
  getOption("tikzDocumentDeclaration"),
  c(
    r"(\renewcommand{\familydefault}{\sfdefault})", 
    r"(\usepackage{helvet})",   # Use sans serif font Helvetica
    r"(\newcommand{\matr}[1]{\mathbf{#1}})",
    r"(\usepackage{xcolor})",
    r"(\definecolor{ggplotRed}{rgb}{0.97255,0.46275,0.42745})",
    r"(\definecolor{ggplotBlue}{rgb}{0,0.74902,0.76863})",
    r"(\definecolor{colorbrewerDarkOrange}{HTML}{d95f02})",
    r"(\definecolor{colorbrewerDarkPurple}{HTML}{7570b3})"
  ))
)
options(tikzLatexPackages = union(
  getOption("tikzLatexPackages"),
  c(
    "\\usepackage{amssymb}",
    "\\usepackage{amsmath}", 
    "\\usepackage{bm}",
    "\\usepackage{graphicx}"
  ))
)

my_get_legend <- function(plot){
  comps <- cowplot::get_plot_component(plot,  "guide-box",return_all = TRUE)
  comps <-  purrr::discard(comps, \(x) ggplot2:::is.zero(x))
  comps[[1]]
}


small_axis <- function(label = NULL, fontsize = 7, arrow_length = 10, label_offset = 1, fix_coord = TRUE, remove_axes = TRUE,
                       arrow_spec = grid::arrow(ends = "both", type = "closed", angle = 20, length = unit(arrow_length / 7, units)),
                       units = "mm", ...){
  coord <- if(fix_coord){
    coord_fixed(clip = "off", ...)
  }else{
    NULL
  }
  axis_theme <- if(remove_axes){
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  }else{
    NULL
  }
  lines <- annotation_custom(grid::polylineGrob(x = unit(c(0, 0, arrow_length), units), y = unit(c(arrow_length, 0, 0), units), 
                                                gp = grid::gpar(fill = "black"),
                                                arrow = arrow_spec))
  text <- if(! is.null(label)){
    annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                     x = unit(label_offset, units), y = unit(label_offset, units), hjust = 0, vjust = 0))
  }else{
    NULL
  }
  list(coord, axis_theme, lines, text)
}

small_arrow <- function(position = c(0.8, 0.95), offset = 0.01, label = NULL, direction = c("x", "y"), 
                        fontsize = 7, arrow_length = unit(10 / 7, "mm"), label_offset = 0, label_hjust = NULL, label_vjust = NULL,
                        arrow_spec = grid::arrow(ends = "last", type = "closed", angle = 20, length = arrow_length),
                        units = "npc"){
  direction <- match.arg(direction)
  if(!grid::is.unit(position)){
    position <- grid::unit(position, units = units)
  }
  if(!grid::is.unit(offset)){
    offset <- grid::unit(offset, units = units)
  }
  if(!grid::is.unit(label_offset)){
    label_offset <- grid::unit(label_offset, units = units)
  }
  if(direction == "x"){
    arrow <- annotation_custom(grid::polylineGrob(x = position, y = c(offset, offset), 
                                                  gp = grid::gpar(fill = "black"),
                                                  arrow = arrow_spec))
    text <- if(! is.null(label)){
      annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                       x = (position[1] + position[2]) / 2, y = offset + label_offset,
                                       hjust = label_hjust, vjust = label_vjust))
    }
  }else{
    arrow <- annotation_custom(grid::polylineGrob(y = position, x = c(offset, offset), 
                                                  gp = grid::gpar(fill = "black"),
                                                  arrow = arrow_spec))
    text <- if(! is.null(label)){
      annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                       y = (position[1] + position[2]) / 2, x = offset + label_offset,
                                       hjust = label_hjust, vjust = label_vjust, rot = 90))
    }
  }
  list(arrow, text)
}

signif_to_zero <- function(x, digits = 6){
  n_signif_digits <- digits - ceiling(log10(abs(x)))
  sign(x) * floor(abs(x) * 10^n_signif_digits) / 10^n_signif_digits
}

######### Custom plotting functions #########

convert_dims <- function( width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1){
  units <- match.arg(units)
  if(units == "inches"){
    units <- "in"
  }
  to_inches <- function(x) x/c(`in` = 1, cm = 2.54, mm = 2.54 * 
                                 10, px = dpi)[units]
  to_inches(c(width, height)) * scale
}

my_pdf <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  grDevices::pdf(filename, width = dim[1], height = dim[2], useDingbats = FALSE, ...)
}


my_tikz <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, stand_alone = TRUE, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  tikzDevice::tikz(filename, width = dim[1], height = dim[2], standAlone = stand_alone, 
                   documentDeclaration = getOption("tikzDocumentDeclaration"), packages = getOption("tikzLatexPackages"), ..., verbose = TRUE)
}

save_plot <- function(filename, plot = ggplot2::last_plot(), width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), dpi = 300, scale = 1, latex_support = FALSE, ...){
  
  old_dev <- grDevices::dev.cur()
  if(latex_support){
    filename <- if(stringr::str_ends(filename, "\\.pdf")){
      paste0(stringr::str_sub(filename, end  = -5L), ".tex")
    }
    my_tikz(filename, width = width, height = height, units = units, dpi = dpi, scale = scale, stand_alone = TRUE)
  }else{
    dim <- convert_dims(width, height, units, dpi, scale)
    dev <- ggplot2:::plot_dev(NULL, filename, dpi = dpi)
    dev(filename = filename, width = dim[1], height = dim[2], ...)
    on.exit(utils::capture.output({
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
    }))
  }
  
  grid::grid.draw(plot)
  
  if(latex_support){
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
    
    withr::with_dir(dirname(filename), {
      tools::texi2pdf(basename(filename), clean = TRUE)
      # Remove .tex file
      file.remove(basename(filename)) 
      raw_file_name <- tools::file_path_sans_ext(basename(filename))
      # rastered images
      ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
      if(length(ras_files) > 0){
        file.remove(ras_files)
      }
    })
  }
  
  invisible(filename)  
}

plot_assemble <- function(..., .plot_objs = NULL, width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), 
                          latex_support = FALSE, show_grid_lines = TRUE, filename = NULL){
  units <- match.arg(units)
  
  plots <- if(is.null(.plot_objs)){
    list(...)
  }else{
    .plot_objs
  }
  
  if(show_grid_lines){
    x_breaks <- scales::breaks_pretty(n = 10)(seq(0, width, length.out = 100))
    y_breaks <- scales::breaks_pretty(n = 10)(seq(0, height, length.out = 100))
  }else{
    x_breaks <- c(0,Inf)
    y_breaks <- c(0,Inf)
  }
  
  if(! is.null(filename)){
    old_dev <- grDevices::dev.cur()
    if(latex_support){
      filename <- if(stringr::str_ends(filename, "\\.pdf")){
        paste0(stringr::str_sub(filename, end  = -5L), ".tex")
      }
      my_tikz(filename, width = width, height = height, units = units, stand_alone = TRUE)
    }else{
      my_pdf(filename, width = width, height = height, units = units)
      on.exit(utils::capture.output({
        grDevices::dev.off()
        if (old_dev > 1) grDevices::dev.set(old_dev)
      }))
    }
  }
  
  
  plotgardener::pageCreate(width = width, height = height, default.units = units, xgrid = diff(x_breaks)[1], ygrid = diff(y_breaks)[1], showGuides = show_grid_lines)
  for(obj in plots){
    if(is.ggplot(obj)){
      plotgardener::plotGG(obj, x = 0, y = 0, width = width, height = height, default.units = units)
    }else if(inherits(obj, "tikz_graphic_fun")){
      grid::grid.draw(obj$FUN(height))
    }else if(grid::is.grob(obj)){
      grid::grid.draw(obj)
    }else if(is.list(obj)){
      stopifnot(! is.null(names(obj)))
      stopifnot("plot" %in% names(obj))
      .x <- obj$x %||% 0
      .y <- obj$y %||% 0
      .width <- obj$width %||% width
      .height <- obj$height %||% height
      .units <- obj$units %||% units
      plotgardener::plotGG(obj$plot, x = .x, y = .y, width = .width, height = .height, default.units = .units)
    }else{
      warning("Cannot handle object of class: ", toString(class(obj)))
    }
  }
  
  if(! is.null(filename)){
    if(latex_support){
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
      withr::with_dir(dirname(filename), {
        tools::texi2pdf(basename(filename), clean = TRUE)
        # Remove .tex file
        file.remove(basename(filename)) 
        raw_file_name <- tools::file_path_sans_ext(basename(filename))
        # rastered images
        ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
        if(length(ras_files) > 0){
          file.remove(ras_files)
        }
      })
    }
  }
}

add_plot <- function(plot, x = 0, y = 0, width = NULL, height = NULL, units = NULL){
  list(plot = plot, x = x, y = y, width = width, height = height, units = units)
}

add_text <- function(label, x = 0, y = 0, fontsize = 12, hjust = 0, vjust = 0, ...){
  list(plot = cowplot::ggdraw() + cowplot::draw_label(label, size = fontsize, hjust = hjust, vjust = vjust, ...), x = x, y = y, width =0, height = 0)
}

# Note that x and y are from the lower left corner (instead of upper left :/)
add_graphic <- function(filename, x = 0, y = 0, width = NULL, height = NULL,
                        units = c("inches", "cm", "mm", "px", "user"),
                        anchor = c("north west", "south west", "base")){
  stopifnot(file.exists(filename))
  units <- match.arg(units)
  anchor <- anchor[1]
  abs_filepath <- tools::file_path_as_absolute(filename)
  size_spec <- if(!is.null(height) && !is.null(width)){
    paste0("[width=", width, units, ", height=", height, units, "]")
  }else if(!is.null(height)){
    paste0("[height=", height, units, "]")
  }else if(! is.null(width)){
    paste0("[width=", width, units, "]")
  }else{
      ""
  }
  content <- paste0(r"(\includegraphics)", size_spec, r"({")", abs_filepath, r"("})")
  res <- list(FUN = (\(figure_height){
    tikzDevice::grid.tikzNode(x = x, y = figure_height - y, units = units, 
                              opts = paste0("draw=none,fill=none,anchor=", anchor),  
                              content = content, draw = FALSE)
  }))
  class(res) <- "tikz_graphic_fun"
  res
}


unnest_named_lists <- function(data, cols, names_to = "name"){
  cols <- tidyselect::eval_select(expr = enquo(cols), data = data,  allow_rename = FALSE)
  map2(cols, names(cols), \(co, na){ 
    val <- data[[co]]
    if(! is.list(val) || is.null(names(val))){
      stop("All columns must be named list columns. Column: ", na , " is not.")
    }
    tibble(..id.. = rep(seq_len(nrow(data)), lengths(val)), ..names.. = unlist(map(val, names)), {{na}} := unlist(val))
  }) |>
    purrr::reduce(dplyr::full_join, by = c("..id..", "..names..")) |>
    (\(x)
     left_join({
       data |>
         dplyr::select(-all_of(cols)) |>
         mutate( ..id.. = row_number())
     }, x, by = "..id.."))() |>
    dplyr::select(-..id.., {{names_to}} := ..names..)
}
