d = evoverse::Clusters(x)

s = 'Set6_42'

sd = d[1:10, ] %>%
  select(id,
         gene,
         # cosmic,
         # function.,
         # mutlocation,
         region,
         vartype,
         starts_with(!!s)) %>%
  separate(id, into = c('chr', 'from', 'to', 'ref', 'alt'), sep = ':')

colnames(sd) = gsub(
  colnames(sd),
  pattern = paste0(s, '[.]'),
  replacement = ''
)

colnames(sd)[ncol(sd):(ncol(sd) - 1)] = c("MOBSTER_VIBER", "MOBSTER")
colnames(sd)[12] = c("VAF_N")

sd$VAF_N = as.numeric(sd$VAF_N)

require(formattable)
require(kableExtra)

sd = sd %>% mutate(
  chr = color_tile("white", "orange")(paste(chr)),
  # VAF = color_bar("indianred3")(VAF)
)

sd = sd %>% column_spec(sd, 6, width = "20em", bold = TRUE, italic = TRUE)


sd  %>%
  # column_spec(1, bold = T) %>%
  # mutate_if(is.numeric, round, digit = 2) %>%
  kable(format = 'html', escape = F, digits = 2) %>%
  add_header_above(c(
    " " = 8,
    "Data" = 4,
    "Clusters" = 2
  )) %>%
  kable_styling(full_width = F) %>%
  pack_rows("Group 1", 4, 7) %>%
  pack_rows("Group 2", 8, 10)


  q# scroll_box(width = "500px", height = "600px") %>%
  save_kable(file = "table1.html", self_contained = T)


  sd[1:10, ] %>%
    # mutate_if(is.numeric, function(x) {
    #   cell_spec(x, bold = T,
    #             color = spec_color(x, end = 0.9),
    #             font_size = spec_font_size(x))
    # }) %>%
    mutate(VAF = cell_spec(
      VAF, color = "white", bold = T,
      background = spec_color(1:10, end = 0.9, option = "A", direction = -1)
    )) %>%
    kable(escape = F, align = "c") %>%
    kable_styling(c("striped", "condensed"), full_width = F)
