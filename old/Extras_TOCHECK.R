#' Title
#'
#' @param data
#' @param cluster
#' @param samples
#' @param palette
#' @param cex
#' @param title
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_mplot_mapping2MOBSTER_clusters = function(data,
                                                  MOBSTER.fits,
                                                  cluster = 'sciClone.cluster',
                                                  samples,
                                                  palette = 'Set2',
                                                  cex = 1,
                                                  tail.color = 'gray',
                                                  title = paste0('Mapping of ', cluster, ' clusters to MOBSTER ones')
)
{
  labels = data[, cluster]
  groups = split(data, f = labels)

  values = lapply(
    groups,
    function(g) {
      lapply(samples, function(s) {

        w = which(g[, paste0('VAF.', s)] > 0)

        t = table(g[w, paste0('MOBSTER.', s, '.cluster')])
        t = as.data.frame(t)

        if(nrow(t)>0) t$Sample = s
        else t = NULL
        t
      })
    })

  values = lapply(values, Reduce, f = rbind)
  values = lapply(seq(values), function(g) {
    cbind(values[[g]], Cluster = names(values)[g])
  })

  values = Reduce(rbind, values)

  # get Palette
  maxBeta = max(sapply(MOBSTER.fits, function(w) w$Kbeta))

  col = scols(
    paste0('C', 1:maxBeta), palette = palette)
  col = c(col, `Tail` = tail.color)


  p1 = ggplot(values,
              aes(Freq, x = Sample, fill = factor(Var1))) +
    geom_bar(stat = "identity",  position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = col, labels = names(col)) +
    theme_light(base_size = 8 * cex) +
    guides(fill=guide_legend(title="MOBSTER Cluster"))  +
    # facet_grid(~Cluster, scales = "free_y", space = "free_y") +
    facet_wrap(~Cluster, scales="free_x", nrow = 1)+
    labs(
      title = bquote(bold(.(title))),
      y = bquote("Number of observations"),
      subtitle = 'Multivariate clustering of read counts'
    ) +
    # geom_hline(yintercept = cutoff, colour = 'red', linetype = "longdash") +
    theme(legend.position="bottom",
          legend.text = element_text(size = 8 * cex),
          legend.key.size = unit(.3 * cex, "cm")
    ) +
    coord_flip()

  dummy <- ggplot(values,
                  aes(Freq, x = Sample)) +
    facet_wrap(~Cluster, scales="free_x", nrow = 1) +
    geom_rect(aes(fill=factor(Cluster)), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_brewer(palette = 'Set2')
  theme_minimal()

  library(gtable)

  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(dummy)

  gtable_select <- function (x, ...)
  {
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
  }

  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip_t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1

  require(grid)
  new_strips <- gtable_select(g2, panels | strips)
  # grid.newpage()
  # grid.draw(new_strips)

  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }

  ## ideally you'd remove the old strips, for now they're just covered
  new_plot <- gtable_stack(g1, new_strips)
  ggpubr::as_ggplot(new_plot)
}

show_VAF_per_segment =  function(segments, muts, cutoff, ...)
{
  plots = lapply(1:nrow(segments),
                 function(s) {
                   # map SNVs to this segment
                   map = map_SNV_2_CNA_segment(segments = segments,
                                               segment.id = s,
                                               muts = muts,
                                               ...)

                   map = map[, startsWith(colnames(map), 'VAF')]

                   map = reshape2::melt(map)
                   map = map[map$value > cutoff,]

                   mxlim = max(1, max(map$value, na.rm = T))

                   # lbl = paste0(
                   #
                   # )

                   ggplot(map, aes(value, fill = variable)) +
                     geom_histogram(binwidth = 0.01) +
                     facet_wrap( ~ variable, scale = 'free', nrow = 1) +
                     coord_cartesian(xlim = c(0, mxlim)) +
                     guides(fill = FALSE) +
                     labs(title = paste(segments[s, ], collapse = ',')) +
                     geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3)
                 })

  figure = ggpubr::ggarrange(
    plotlist = plots,
    nrow = nrow(segments),
    ncol = 1
  )

  figure
}

#' Title
#'
#' @param fit_with_MOBSTER
#' @param fit_without_MOBSTER
#' @param CCF.cutoff
#' @param title
#' @param cex_legend
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_clustering_assignment_summary = function(fit_with_MOBSTER,
                                                     fit_without_MOBSTER,
                                                     CCF.cutoff = .70,
                                                     title = fit_with_MOBSTER$description,
                                                     cex_legend = 1)
{
  x = fit_with_MOBSTER
  y = fit_without_MOBSTER

  pio::pioTit("Plotting ...")

  pio::pioStr("Cutoff to binarize adjusted VAF",
              paste0(CCF.cutoff / 2, " (CCF ", CCF.cutoff, ')'))

  CCF.cutoff = CCF.cutoff / 2

  sample = x$samples[1]
  k = length(x$samples)

  # mutations = full_join(VAF_table(x), MOBSTER_clusters, by = 'id')

  # M_clusters = MClusters(x)
  # colnames(M_clusters)[1+1:k] = paste0("MOBSTER.", colnames(M_clusters)[1+1:k])
  #
  # mutations = full_join(VAF_table(x), M_clusters, by = 'id')
  #
  # mutations = VAF(x) %>% select(-id)
  # colnames(mutations) = sub('.VAF',  '', colnames(mutations))
  #
  # pheatmap::pheatmap(mutations)
  #
  # bin.mutations = as.data.frame( %>% select(-id))
  # rownames(bin.mutations) = VAF_table(x)$id
  # bin.mutations[bin.mutations > 0] = 1



  suffix = paste0(" >", CCF.cutoff)

  # VAF data, to which we add the median VAF (for sorting)
  VAF_data = VAF_table(y)
  VAF_data$VAF_median = apply(VAF_data[1 + 1:k], 1, mean)

  # Binarize data from VAF
  is = function(x) {
    is.numeric(x) &
      x > CCF.cutoff
  } # We binarize to 1 everything above CCF.cutoff

  bin.mutations = VAF_data %>%
    select(id, ends_with('.VAF')) %>%
    mutate_if( ~ any(is(.x)),  ~ if_else(is(.x), 1, 0))

  # We count also how many =1 binarized mutations per sample are there
  bin.mutations$Above_cutoff = rowSums(bin.mutations[1 + 1:k])

  # We assemble mutations, with binary data as well
  mutations = full_join(bin.mutations,
                        VAF_data,
                        by = 'id',
                        suffix = c(suffix, ""))

  # Prepare clustering assignments from MOBSTER, and Binomial data
  M_clusters = MClusters(x)
  colnames(M_clusters)[1 + 1:k] = gsub(pattern = 'cluster.',
                                       replacement = 'MOBSTER_',
                                       colnames(M_clusters)[1 + 1:k])

  # paste0(colnames(M_clusters)[1+1:k], '.MOBSTER')
  B_clusters = BClusters(x) %>% rename(Binomial = cluster.Binomial)

  # All clustering outputs
  all_clusters = full_join(B_clusters, M_clusters, by = 'id') %>%
    select(starts_with('Binomial'),
           starts_with('MOBSTER'),
           'id')

  is.tail = apply(all_clusters, 1, function(w)
    any(w == 'Tail', na.rm = TRUE))

  all_clusters$tail_mutation = paste(is.tail)

  mutations = full_join(all_clusters, mutations, by = 'id')
  # mutations = mutations %>% arrange(desc(k))

  cl_cols = colnames(mutations)[startsWith(colnames(mutations), 'MOBSTER')]
  vf_cols = colnames(mutations)[endsWith(colnames(mutations), '.VAF')]

  # View(mutations %>% mutate(medianVAF = median(!!vf_cols)))

  # mutations = mutations %>% arrange(
  #   desc(k),
  #   `Binomial`,
  #   `MOBSTER Set7_55`,
  #   `MOBSTER Set7_57`,
  #   `MOBSTER Set7_59`,
  #   `MOBSTER Set7_62`
  # )

  mutations = mutations %>% arrange_(.dots = c('desc(Above_cutoff)', 'Binomial', cl_cols))

  mutations = mutations %>%
    group_by_(.dots =  c('Above_cutoff', 'Binomial', cl_cols)) %>%
    arrange(VAF_median, .by_group = TRUE) %>%
    ungroup()


  # mutations = mutations[
  #   with(mutations, order(cl_cols)),
  #   ]

  require(pheatmap)

  core_pheat = as.data.frame(mutations %>% select(ends_with('.VAF')))
  rownames(core_pheat) = mutations$id



  core_pheat_colors = rev(mobster:::scols(1:50, 'Spectral'))

  core_pheat_colors = viridisLite::viridis(
    50,
    alpha = 1,
    begin = 0,
    end = 1,
    direction = 1,
    option = "D"
  )
  core_pheat_colors[1] = 'gray95'

  # viridisLite::viridis(50, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
  M = max(mutations %>% select(ends_with('.VAF')))
  bin_color = M / 50
  core_pheat_colors[floor(CCF.cutoff / bin_color)] = 'red'

  M = max(mutations %>% select(ends_with('.VAF')), na.rm = TRUE)
  bin_color = M / 50
  above = floor(abs((M - CCF.cutoff) / 0.02))
  below = floor(abs(CCF.cutoff / 0.02))

  core_pheat_colors =
    c('gray95',
      (mobster:::scols(1:below, 'YlGnBu')),
      'red',
      mobster:::scols(1:above, 'Oranges'))

  col_below = mobster:::scols(1:below, 'YlGnBu')
  core_pheat_colors =
    c('gray95', col_below, rep(col_below[length(col_below)], above))
  core_pheat_colors[floor(CCF.cutoff / (M / 35))] = 'red'




  # Annotate binarized VAF
  annotations_binarized_vaf = as.data.frame(mutations %>% select(ends_with(suffix), Above_cutoff))
  rownames(annotations_binarized_vaf) = mutations$id

  # annotations_binarized_vaf[which(annotations_binarized_vaf == 1, arr.ind = T)] = "TRUE"
  # annotations_binarized_vaf[which(annotations_binarized_vaf != "TRUE", arr.ind = T)] = "FALSE"

  # Annotate clusters
  annotations_clusters = as.data.frame(mutations %>% select(Binomial, starts_with('MOBSTER'), tail_mutation))
  rownames(annotations_clusters) = mutations$id

  # Color Binomial clusters
  colors_B_clusters = mobster:::scols(unique(annotations_clusters$Binomial), 'Set1')

  colors_B_clusters[is.na(names(colors_B_clusters))] = 'gray'

  annotation_colors = list(Binomial = colors_B_clusters)



  # Color MOBSTER clusters
  for (sample in x$samples) {
    # assgn = paste0('cluster.', sample, '.MOBSTER')
    assgn = paste0('MOBSTER_', sample)

    lbl = unique(annotations_clusters[, assgn])
    lbl = sort(lbl, na.last = T)

    pal = rev(wesanderson::wes_palette("FantasticFox1"))
    colors_M_clusters = mobster:::scols(lbl, 'Set2')

    colors_M_clusters = pal[1:length(colors_M_clusters)]
    names(colors_M_clusters) = lbl

    colors_M_clusters[is.na(names(colors_M_clusters))] = 'gray'
    colors_M_clusters[names(colors_M_clusters) == 'Tail'] = 'darkgray'

    l = list(colors_M_clusters)
    names(l) = assgn

    annotation_colors = append(annotation_colors, l)
  }

  tail_yesno_color = list(c(`TRUE` = 'darkgray', `FALSE` = 'white'))
  names(tail_yesno_color) = 'tail_mutation'

  annotation_colors = append(annotation_colors, tail_yesno_color)

  annotations_binarized_vaf = annotations_binarized_vaf[mutations$id, , drop = F]
  annotations_clusters = annotations_clusters[mutations$id, , drop = F]


  gaps_Above_cutoff = unique(annotations_binarized_vaf$Above_cutoff)
  gaps_Above_cutoff = lapply(gaps_Above_cutoff, function(s)
    which(annotations_binarized_vaf$Above_cutoff == s))
  gaps_Above_cutoff = sapply(gaps_Above_cutoff, function(w)
    w[1])

  gaps_Above_cutoff = gaps_Above_cutoff[!is.na(gaps_Above_cutoff)]

  p = pheatmap(
    main = title,
    core_pheat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    color = core_pheat_colors,
    gaps_row = gaps_Above_cutoff,
    # annotation_legend = FALSE,
    annotation_row = cbind(annotations_binarized_vaf, annotations_clusters),
    annotation_colors = annotation_colors,
    cellwidth = 20,
    border_color = NA,
    silent = TRUE
  )

  # library(grid)

  # grid.gedit("GRID.rect.531", gp = gpar(cex=.1))
  # if(cex_legend != 1)
  # grid::grid.gedit("annotation_legend", gp = gpar(cex=cex_legend))
  p$gtable$grobs[[6]] = grid::editGrob(p$gtable$grobs[[6]], gp = gpar(cex =
                                                                        cex_legend))

  p
}

