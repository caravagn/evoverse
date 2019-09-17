#' Plot the VAF per kariotype.
#'
#' @description Plot the adjusted VAF histogram from the mutations
#' in \code{x}, split by kariotype.
#'
#' @param x An object of class \code{mbst_data}.
#' @param range The VAF range (left to right), whichi is by default [0;+inf).
#' @param cex Cex of the plot.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' data("example_mvmobster")
plot_VAF_kariotypes = function(x,
                               range = c(0, Inf),
                               cex = 1)
{
  stopifnot(inherits(x, "mbst_data"))

  values = x$VAF_cn_adjustment %>%
    mutate(Genotype = paste0(Major, ':', minor))

  ggplot(values, aes(adj_VAF, fill = sample)) +
    geom_histogram(binwidth = 0.01) +
    facet_grid(Genotype ~ sample, scales = 'free') +
    theme_light(base_size =  8 * cex) +
    guides(fill = FALSE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    geom_vline(
      aes(xintercept = 0.5),
      colour = 'red',
      linetype = "longdash",
      size = .3
    ) +
    geom_vline(
      aes(xintercept = 0.25),
      colour = 'red',
      linetype = "longdash",
      size = .3
    ) +
    geom_vline(
      aes(xintercept = 0.33),
      colour = 'blue',
      linetype = "longdash",
      size = .3
    ) +
    geom_vline(
      aes(xintercept = 0.66),
      colour = 'blue',
      linetype = "longdash",
      size = .3
    ) +
    geom_vline(
      aes(xintercept = 1),
      colour = 'black',
      linetype = "longdash",
      size = .2
    ) +
    xlim(range[1], range[2]) +
    labs(
      title = bquote(bold("Adjusted VAF per kariotype")),
      subtitle = paste0('Plot range: ', range[1], '-', range[2]),
      x = 'Adjusted VAF',
      y = "Counts"
    )
}

#' Plot the counts of each kariotype.
#'
#' @description Plot the number of segments per sample and
#' per kariotype, in \code{x}, as a barplot.
#'
#' @param x An object of class \code{mbst_data}.
#' @param cex Cex of the plot.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' data("example_mvmobster")
plot_count_kariotypes = function(x, cex = 1)
{
  stopifnot(inherits(x, "mbst_data"))

  values = x$VAF_cn_adjustment %>%
    mutate(Genotype = paste0(Major, ':', minor)) %>%
    select(sample, Genotype) %>%
    group_by(sample, Genotype) %>%
    filter(row_number() == 1) %>%
    summarise(count = n()) %>%
    ungroup()

  ggplot(values, aes(y = count, x = Genotype, fill = sample)) +
    geom_bar(stat = 'identity') +
    facet_wrap( ~ sample, scales = 'free', nrow = 1) +
    theme_light(base_size =  8 * cex) +
    guides(fill = FALSE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    labs(title = bquote(bold("Number of segments per kariotype")),
         x = 'Kariotype',
         y = 'Number of segments')
}

#' @description Plot the number of mutations per sample and
#' per kariotype, in \code{x}, as a barplot.
#'
#' @param x An object of class \code{mbst_data}.
#' @param cex Cex of the plot.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("example_mvmobster")
plot_count_muts_per_kariotype = function(x, cex = 1)
{
  stopifnot(inherits(x, "mbst_data"))

  values = x$VAF_cn_adjustment %>%
    mutate(Genotype = paste0(Major, ':', minor)) %>%
    select(sample, Genotype) %>%
    group_by(sample, Genotype) %>%
    summarise(count = n()) %>%
    ungroup()

  ggplot(values, aes(y = count, x = Genotype, fill = sample)) +
    geom_bar(stat = 'identity') +
    facet_wrap( ~ sample, scales = 'free', nrow = 1) +
    theme_light(base_size = 8 * cex) +
    guides(fill = FALSE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(.3, "cm")
    ) +
    labs(title = bquote(bold("Number of mutations per kariotype")),
         x = 'Kariotype',
         y = 'Number of mutations')
}





# load("/Users/gcaravagna/Downloads/rCGH/data/hg19.rda")
# hg19 = as_tibble(hg19)
# hg19 = hg19 %>% mutate(from = cumlen)
# hg19 = hg19 %>% mutate(to = from + length)
# hg19 = hg19 %>% mutate(centromerStart = from + centromerStart)
# hg19 = hg19 %>% mutate(centromerEnd = from + centromerEnd)
# hg19 = hg19 %>% mutate(chr = paste0('chr', chrom))
# hg19 = hg19 %>% select(chr, length, from, to, centromerStart, centromerEnd)
# hg19$chr[hg19$chr == 'chr23'] = 'chrX'
# hg19$chr[hg19$chr == 'chr24'] = 'chrY'
# chr_coordinate_hg19 = hg19
# save(chr_coordinate_hg19, file = 'data/chr_coordinate_hg19.RData')


#' Plot copy number segments.
#'
#' @description
#'
#' Plots the copy number segments available for a required
#' set of samples. Returns a list of plots named after each
#' input sample via \code{samples}. For every sample it plots
#' the minor and major segments, highlights diploid regions
#' (kariotype minor = 1, Major = 1) and annotates on top the
#' number, the VAF and the DP of each annotated mutation. Each
#' returned object is a figure arranged via \code{ggarrange}.
#'
#' @param x An object of class \code{mbst_data}.
#' @param samples Required samples.
#' @param ... Extra parameters passed to the internal
#' hidden function \code{plt_CNsegments}.
#'
#' @return A list of \code{ggarrange} objects named after samples.
#'
#' @export
#'
#' @examples
#' TODO
plot_CNsegments = function(x,
                           samples = x$samples,
                           ...)
{
  pl = lapply(samples,
              function(w)
                plt_CNsegments(x, sample = w, ...))
  names(pl) = samples

  pl
}


plt_CNsegments = function(x, sample, chromosomes = paste0('chr', c(1:22, 'X', 'Y')))
{
  data('chr_coordinates_hg19', package = 'mobster')

  chr_coordinates_hg19 = chr_coordinates_hg19 %>% filter(chr %in% chromosomes)

  low = min(chr_coordinates_hg19$from)
  upp = max(chr_coordinates_hg19$to)

  stop("Adjust Muts for absolute coordinates.")

  baseplot = ggplot(chr_coordinates_hg19) +
    theme_classic() +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .3,
      colour = 'gainsboro'
    ) +
    geom_vline(xintercept = chr_coordinates_hg19$from,
               size = 0.3,
               colour = 'black') +
    geom_label(
      data = chr_coordinates_hg19,
      aes(
        x = chr_coordinates_hg19$from,
        y = -0.5,
        label = gsub('chr', '', chr_coordinates_hg19$chr)
      ),
      hjust = 0,
      colour = 'white',
      fill = 'black',
      size = 3
    ) +
    geom_hline(yintercept = 0,
               size = 1,
               colour = 'gainsboro') +
    geom_hline(
      yintercept = 1,
      size = .3,
      colour = 'black',
      linetype = 'dashed'
    ) +
    labs(x = "Coordinate",
         y = "Major/ minor allele",
         caption = x$description) +
    ggpubr::rotate_y_text() +
    xlim(low, upp)

  # Segments regions
  segments = x$segments %>%
    filter(sample == !!sample, chr %in% chromosomes) %>%
    spread(variable, value) %>% arrange(chr)

  # if there are 1-1 segments, shadow them
  one_one = segments %>% filter(Major == 1, minor == 1)
  if (nrow(one_one) > 0)
    baseplot = baseplot +
    geom_rect(
      data = one_one,
      aes(
        xmin = from,
        xmax = to,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = .2,
      fill = 'forestgreen'
    )

  # Segments
  baseplot = baseplot +
    geom_segment(
      data = segments,
      aes(
        x = from,
        xend = to,
        y = Major,
        yend = Major
      ),
      size = 1.5,
      colour = 'red'
    ) +
    geom_segment(
      data = segments %>% mutate(minor = minor - 0.1),
      aes(
        x = from,
        xend = to,
        y = minor,
        yend = minor
      ),
      size = 1,
      colour = 'steelblue'
    )

  # Scale size
  if (max(segments$Major) < 5)
    baseplot = baseplot + ylim(-0.5, 5)

  # Mutation VAF
  mutations = x$map_mut_seg %>%
    filter(chr %in% chromosomes)
  vaf = VAF(x, ids = mutations %>% pull(id), samples = sample) %>% filter(value > 0) %>% pull(id)
  mutations = mutations %>%
    filter(id %in% vaf)

  vaf_mutations = VAF(x, ids = mutations %>% pull(id), samples = sample)
  vaf_mutations = vaf_mutations %>% left_join(mutations, by = 'id')

  quant = quantile(vaf_mutations$value, probs = c(.1, .99))
  vaf_mutations = vaf_mutations %>% filter(value > quant[1], value < quant[2])

  dpplot = ggplot(vaf_mutations, aes(x = from, y = value, color = value)) +
    stat_density_2d(
      aes(fill = ..density..),
      geom = "raster",
      contour = FALSE,
      alpha = 1,
      n = 100
    ) +
    scale_fill_distiller(palette = "Spectral") +
    theme_classic() +
    xlim(low, upp) +
    # ylim()
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'Den.') +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 1), limits = c(0, 1)) +
    # ylim(0, 1) +
    guides(fill = FALSE)

  pplot = ggplot(vaf_mutations, aes(x = from, y = value, color = value)) +
    geom_point(size = 1e-2, alpha = .3) +
    # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha = 1) +
    # scale_fill_distiller(palette= "Spectral") +
    geom_hline(
      yintercept = 0.5,
      col = 'forestgreen',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = 0.25,
      col = 'steelblue',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = 0.5 + 0.25,
      col = 'darkred',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = 0,
      col = 'black',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = 1,
      col = 'black',
      size = .3,
      linetype = 'dashed'
    ) +
    theme_classic() +
    xlim(low, upp) +
    # scale_x_continuous(breaks = c(low, upp)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'VAF') +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 1),
                       limits = c(0, 1)) +
    # ylim(0, 1) +
    scale_color_gradientn(
      colours = c('steelblue', 'forestgreen', 'darkred'),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    guides(colour = FALSE)

  # Mutation coverage
  mutations = x$map_mut_seg %>%
    filter(chr %in% chromosomes)
  depth = DP(x, ids = mutations %>% pull(id), samples = sample) %>% filter(value > 0) %>% pull(id)
  mutations = mutations %>%
    filter(id %in% depth)

  depth_mutations = DP(x, ids = mutations %>% pull(id), samples = sample)
  depth_mutations = depth_mutations %>% left_join(mutations, by = 'id')

  q_cov = quantile(depth_mutations$value, probs = c(.01, .99))
  names(q_cov) = q_cov

  pplot_depth = ggplot(depth_mutations, aes(x = from, y = value, color = value)) +
    geom_point(size = 1e-2, alpha = .3) +
    geom_hline(
      yintercept = median(depth_mutations$value),
      col = 'black',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = q_cov[1],
      col = 'black',
      size = .3,
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = q_cov[2],
      col = 'black',
      size = .3,
      linetype = 'dashed'
    ) +
    # geom_hline(yintercept = 0.25, col = 'steelblue', size = .3, linetype = 'dashed') +
    # geom_hline(yintercept = 0.5+0.25, col = 'darkred', size = .3, linetype = 'dashed') +
    theme_classic() +
    xlim(low, upp) +
    # scale_x_continuous(breaks = c(low, upp)) +
    scale_y_continuous(breaks = q_cov, limits = q_cov) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'DP') +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_color_gradientn(
      colours = c('steelblue', 'orange', 'darkred'),
      breaks = c(q_cov[1], median(depth_mutations$value), q_cov[2]),
      limits = q_cov
    ) +
    guides(colour = FALSE)

  # Histogram of mutation counts with 1 megabase bins
  binsize = 1e6

  hplot = ggplot(mutations, aes(x = from)) +
    geom_histogram(aes(y = ..count..), binwidth = binsize, fill = 'black')

  m = max(ggplot_build(hplot)$data[[1]]$count) + 20
  hplot = hplot +
    theme_classic() +
    xlim(low, upp) +
    # scale_x_continuous(low, upp - upp*binsize) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggpubr::rotate_y_text() +
    # scale_y_continuous(
    #   breaks = c(0, m+1, m + 10), labels = c(0, m, '')) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2),
                       limits = c(0, m)) +
    labs(y = 'n')
  # geom_hline(yintercept = m, linetype = 'dashed', size = 0.3, color = 'gray') +
  # annotate("text", x = 0, y = m*0.9, label = ceiling(m), size = 2)

  figure = ggarrange(
    hplot,
    pplot,
    pplot_depth,
    baseplot,
    nrow = 4,
    ncol = 1,
    heights = c(.2, .2, .2, 1)
  )

  annotate_figure(figure,
                  top = text_grob(
                    bquote(
                      bold('   ' ~ .(sample)) ~
                        ' -  Copy Number profile (' * .(nrow(segments)) ~ 'seg., n =' ~ .(nrow(mutations)) * ')'
                    ),
                    color = "black",
                    face = "bold",
                    hjust = 0,
                    x = 0
                  ))
}


#' Preview fiters effect on mutations
#'
#' @description mvMOBSTER has functions to subset the mutations
#' that one has annotated in the data. With this function one
#' can generate plots that show which mutations would be filtered
#' for a certain set of parameters of the filters. This allows
#' one to tune the filters available in mvMOBSTER.
#'
#' @param x An object of class \code{mbst_data}.
#' @param VAF_min VAF values below this parameter are considered as 0 (undetected);
#' see function \code{\link{mobster_flt_minvaf}}. This filter applies to each
#' sample independently.
#' @param NV_min NV values below this parameter are considered as 0 (undetected);
#' see function \code{\link{mobster_flt_minnv}}. This filter applies to each
#' sample independently.
#' @param min.DP The depth value for each mutation must be higher than \code{min.DP},
#' whenever it is greater than 0 in a sample; see function \code{\link{mobster_flt_dprange}}.
#' @param max.DP The depth value for each mutation must be lower than \code{max.DP};
#' see function \code{\link{mobster_flt_dprange}}.
#' @param min.VAF The VAF value for each mutation must be higher than \code{min.VAF},
#' whenever it is greater than 0 in a sample; see function \code{\link{mobster_flt_vafrange}}.
#' @param max.VAF The depth value for each mutation must be lower than \code{min.VAF};
#' see function \code{\link{mobster_flt_vafrange}}.
#' @param x.lim X-axis limits for the plot (e.g., \code{c(0,1)}); NA is no limits.
#' @param y.lim Y-axis limits for the plot (e.g., \code{c(0,1)}); NA is no limits.
#'
#' @return A list of \code{ggpubr} figures (per filter, and per sample) with all annotated mutations
#' coloured according to their status (if they would be removed by a filter, or not).
#'
#' @export
#'
#' @examples
#' TODO
plot_filters = function(x,
                        VAF_min,
                        NV_min,
                        min.DP,
                        max.DP,
                        min.VAF,
                        max.VAF,
                        x.lim = NA,
                        y.lim = NA)
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Pairwise plotting function
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-
  aux_fun_plt_filters = function(obj, x, y)
  {
    pp = VAF_table(obj, samples = c(x, y), suffix = '')

    #
    # FILTER #1 - VAF > minimum value.
    #
    pp$filter_minVAF =
      apply(pp %>% select(-id),
            1,
            function(w)
            {
              w1 = (w[1] > 0) & (w[1] < VAF_min)
              w2 = (w[2] > 0) & (w[2] < VAF_min)

              if (!w1 & !w2)
                return('No')
              else
                return("Yes")
            })

    #
    # FILTER #2 - Depth across min and max in all samples
    #
    ids = DP(obj) %>%  filter(value > 0 &
                                (value < min.DP |
                                   value > max.DP)) %>%
      select(id) %>%
      distinct() %>% pull(id)

    pp = pp %>%
      mutate(filter_dprange = ifelse(id %in% !!ids, 'Yes', 'No'))

    #
    # FILTER #3 - NV > minimum value.
    #
    ppNV = NV_table(obj, samples = c(x, y), suffix = '')

    ppNV$filter_minNV =
      apply(ppNV %>% select(-id),
            1,
            function(w)
            {
              w1 = (w[1] > 0) & (w[1] < NV_min)
              w2 = (w[2] > 0) & (w[2] < NV_min)

              if (!w1 & !w2)
                return('No')
              else
                return("Yes")
            })

    pp$filter_minNV = "No"
    pp[ppNV %>% filter(filter_minNV == 'Yes') %>% pull(id), 'filter_minNV'] = 'Yes'


    #
    # FILTER #4 - NV > minimum value.
    #
    ids = VAF(obj) %>%  filter(value > 0 &
                                 (value < min.VAF |
                                    value > max.VAF)) %>%
      select(id) %>%
      distinct() %>% pull(id)

    pp = pp %>%
      mutate(filter_VAFrange = ifelse(id %in% !!ids, 'Yes', "No"))

    pl = function(attr, t, cex = 1, sub) {
      g = ggplot(pp,
                 aes(
                   x = eval(parse(text = x)),
                   y = eval(parse(text = y)),
                   color = eval(parse(text = attr))
                 )) +
        geom_point(size = .5 * cex) +
        theme_light() +
        labs(
          x = x,
          y = y,
          title = t,
          subtitle = sub
        ) +
        guides(color = guide_legend("Filter")) +
        scale_color_manual(values = c(
          `Yes` = 'red',
          `No` = 'black',
          `NA` = 'gray'
        ))

      if (!is.na(x.lim))
        g = g + xlim(x.lim[1], x.lim[2])
      if (!is.na(y.lim))
        g = g + ylim(y.lim[1], y.lim[2])

      g
    }

    p1 = pl(attr = 'filter_minVAF',
            paste0('VAF > ', VAF_min),
            sub = "Filter type: per sample.")
    p2 = pl(attr = 'filter_dprange',
            paste0('DP range [', min.DP, '; ', max.DP , ']'),
            sub = "Filter type: all samples.")
    p3 = pl(attr = 'filter_minNV',
            paste0('NV > ', NV_min),
            sub = "Filter type: per sample.")
    p4 = pl(attr = 'filter_VAFrange',
            paste0('VAF range [', min.VAF, '; ', max.VAF , ']'),
            sub = "Filter type: all samples")

    pp = pp %>% mutate(
      flt = (filter_minVAF == 'Yes') |
        (filter_dprange == 'Yes') |
        (filter_minNV == 'Yes') | (filter_VAFrange == 'Yes'),
      filter = ifelse(flt, "Yes", "No")
    )

    pall = pl(attr = 'filter', paste0('Filter (any)'),  sub = "Mutations triggering at least one filter.") + facet_wrap( ~
                                                                                                                           filter)
    # Arrange the plots
    ggarrange(
      p1,
      p3,
      p2,
      p4,
      pall,
      nrow = 1,
      ncol = 5,
      common.legend = TRUE,
      legend = 'right'
    )
  }

  pairs = combn(x$samples, 2)

  pl = apply(pairs,
             2,
             function(w) {
               aux_fun_plt_filters(obj = x, x = w[1], y = w[2])
             })
  names(pl) = apply(pairs, 2, paste, collapse = '_vs_')

  pl
}

