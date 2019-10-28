data("example_evoverse")

x = mobster::Clusters(example_evoverse$fit_MOBSTER$Set7_55$best) %>%
  filter(row_number() < 100)

require(kableExtra)

x %>%
  kable() %>%
  kable_styling()


smp_col = function(w, palette)
{
  labels = w %>% unique
  nlabels = labels %>% length

  colors = RColorBrewer::brewer.pal(n = nlabels, name = palette)
  names(colors) = labels

  colors[labels]
}

cluster_colors = smp_col(x %>% pull(cluster), 'Set1')
kario_colors = smp_col(x %>% pull(karyotype), 'Accent')
sample_colors = smp_col(x %>% pull(sample), 'Set2')


x[1:20, ]  %>%
  mutate(
    VAF = cell_spec(
      VAF,
      color = "white",
      bold = T,
      background = spec_color(
        1:10,
        end = 0.9,
        option = "A",
        direction = -1
      )
    ),
    karyotype = cell_spec(karyotype, color = kario_colors[karyotype]),
    sample = cell_spec(sample_colors, color = sample_colors[sample]),
    cluster = cell_spec(cluster, color = cluster_colors[cluster]),
    Tail = formattable::color_bar("black")(round(Tail, 1)),
    C1 = formattable::color_bar("black")(round(C1, 1))
  ) %>%
  kable(escape = F, align = "l") %>%
  kable_styling(c("striped", "condensed"), full_width = T)
