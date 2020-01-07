

x = Set7

calls = lapply(x$samples,
               function(s)
               {
                 x$CNAqc[[s]]$cna %>%
                 mutate(
                   label = paste(Major, minor, sep = ':'),
                   CN = minor + Major,
                   sample = s
                   ) %>%
                 select(chr, from, to, label, CN, sample)
                 })

calls = lapply(calls, CNAqc:::relative_to_absolute_coordinates)

calls_flat = Reduce(full_join, calls) %>%
  mutate(
    label = ifelse(label %in% names(KARYO_colors), label, 'other')
  )

chromosomes = calls_flat$chr %>% unique

data("chr_coordinates_hg19", package = "CNAqc")
chr_coordinates_hg19 = chr_coordinates_hg19 %>% filter(chr %in% chromosomes)


low = min(chr_coordinates_hg19$from)
upp = max(chr_coordinates_hg19$to)

bl_genome =
  ggplot() +
  CNAqc:::my_ggplot_theme() +
  geom_vline(xintercept = chr_coordinates_hg19$centromerStart,
             size = 0.1,
             colour = "gainsboro") +
  geom_vline(xintercept = chr_coordinates_hg19$centromerEnd,
             size = 0.1,
             colour = "gainsboro") +
  geom_vline(xintercept = chr_coordinates_hg19$from,
             size = 0.2,
             colour = "black") +
  labs(x = "Location", y = "Sample") +
  scale_x_continuous(
    breaks = c(0, chr_coordinates_hg19$from, upp),
    labels = c("", gsub(pattern = 'chr', replacement = '', chr_coordinates_hg19$chr), "")
    )

KARYO_colors =
  c(
    `0:0` = 'darkblue',
    `1:0` = 'steelblue',
    `1:1` = 'palegreen3',
    `2:0` = 'turquoise4',
    `2:1` = 'orange',
    `2:2` = 'indianred3',
    `other` = 'darkred'
  )


seg_id = pio:::nmfy(x$samples, seq_along(x$samples))
calls_flat$sample_id = seg_id[calls_flat$sample]

calls %>% mutate()

bl_genome +
  geom_segment(
    data = calls_flat,
    aes(
      x = from,
      xend = to,
      y = sample_id,
      yend = sample_id,
      color = label
    ),
    size = 5
  ) +
  scale_color_manual(values = KARYO_colors) +
  coord_polar(theta = 'x', clip = 'off') +
  guides(color = guide_legend('Karyotype')) +
  theme(legend.key.height = unit(.1, "cm")) +
  labs(title = x$description) +
  ylim(0, max(seg_id) + 1)




