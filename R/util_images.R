# Annotate title on top of a plot
img_annot_t = function(x, t, ...) {
  ggpubr::annotate_figure(x,
                          top = ggpubr::text_grob(
                            t,
                            hjust = 0,
                            x = 0,
                            ...
                          ))
}


# Add margins to a plot
img_add_margin = function(x, fact, w = 1, h = 1) {
  x + theme(plot.margin = unit(c(1 * h, 1 * w, 1 * h, 1 * w) * fact, "cm"))
}


# Button like
button = function(color, radio, size, font = 'plain')
{
  ggplot(data = data.frame(
    x = 0,
    y = 0,
    label = radio,
    stringsAsFactors = FALSE
  ),
  aes(x = x , y = y, label = label)) +
    theme_void() +
    theme(plot.background = element_rect(fill = color)) +
    geom_text(fontface = font, hjust = 0)
}
