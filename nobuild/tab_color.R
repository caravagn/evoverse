data("example_evoverse")

# load('~/Documents/GitHub/test.dbpmm/Real Data/[Used] Koerber_et_al/H043-B7R7/post_dpl_fit.RData')

x = mobster::Clusters(example_evoverse$fit_MOBSTER$Set7_55$best) %>%
  sample_n(100)

# x = mobster::Clusters(post_dpl_fit$fit$best) %>%
#   sample_n(100) %>%
#   rename

require(kableExtra)

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

library(knitr)
library(kableExtra)
library(formattable)
library(dplyr)

x  %>%
  arrange(desc(C1)) %>%
  mutate(
    VAF = case_when(
      (VAF < 0.2) ~ cell_spec(VAF, "html", color = "darkblue"),
      (VAF >= 0.2 & VAF < 0.4) ~ cell_spec(VAF, "html", color = "steelblue"),
      (VAF >= 0.4 & VAF < 0.6) ~ cell_spec(VAF, "html", color = "orange"),
      (VAF >= 0.6 & VAF < 0.8) ~ cell_spec(VAF, "html", color = "red"),
      (VAF > 0.8) ~ cell_spec(VAF, "html", color = "darkred"),
      TRUE ~ cell_spec(VAF, "html", color = "black")
    ),
    # VAF = cell_spec(
    #   round(VAF, 2),
    #   color = "white",
    #   bold = T,
    #   background = spec_color(
    #     1:10,
    #     end = 0.9,
    #     # option = "B",
    #     direction = -1
    #   )
    # ),
    karyotype = cell_spec(karyotype, bold = T, color = 'white', background = kario_colors[karyotype]),
    # karyotype = color_tile('white', kario_colors[karyotype])(karyotype),
    sample = cell_spec(sample, bold = T, color = 'white', background = sample_colors[sample]),
    cluster = cell_spec(cluster, bold = T, color = 'white', background = cluster_colors[cluster]),
    Tail = formattable::color_bar(cluster_colors['Tail'])(round(Tail, 1)),
    C1 = formattable::color_bar(cluster_colors['C1'])(round(C1, 1))
  ) %>%
  kable(escape = F, align = "l") %>%
  kable_styling(c("striped", "condensed"), full_width = T)
