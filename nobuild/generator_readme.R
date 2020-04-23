PKG = c("CNAqc", 'ctree', 'mtree', 'BMix', 'VIBER')

DPKG = c(
  `CNAqc` = 'Copy Number Alterations quality check package to assess the concordance between segments and somatic mutations, peak detection metrics and Cancer Cell Fractions estimation.',
  `ctree` = 'Cancer clones trees from multi-region bulk sequencing data, built of Cancer Cell Franctions (CCFs) clusters computed by tumour subclonal deconvolution.',
  `mtree` = 'Cancer mutation trees from multi-region bulk or single-cell sequencing data, where  the presence or absence of a somatic mutation, a Copy Number event, or any other event is available in binary format.',
  `BMix` = 'Univariate Binomial and/ or Beta-Binomial mixture models for clustering of sequencing read counts. To be used for subclonal deconvolution from bulk sequencing.',
  `VIBER` = 'Multivariate Binomial  mixture models fit by variational inference, which can be ued to cluster read counts from multi-region sequencing data and carry out subclonal deconvolution from bulk sequencing.'
  )

TPKG = c(
  `CNAqc` = 'Copy Number Alterations quality check',
  `ctree` = 'Cancer clones trees',
  `mtree` = 'Cancer mutation trees',
  `BMix` = 'Univariate Binomial and Beta-Binomial mixture models',
  `VIBER` = 'Variational multivariate Binomial mixture models'
)

# img src
isrc = function(x)
{
  paste0('<img src=\"https://caravagn.github.io/', x, '/reference/figures/logo.png\" width="8%">')
}

# img src align=right
isrcr = function(x)
{
  paste0('<img src=\"https://caravagn.github.io/', x, '/reference/figures/logo.png\" align="right" width="8%">')
}

# A href
ahr = function(x, y)
{
  paste0('<a href=\"https://caravagn.github.io/', x, '\" width="8%">', y, '</a>', collapse = '')
}

# stripe
generator_stripe = function(PKG)
{
  block = sapply(
    PKG,
    function(x){
      ahr(x, isrc(x))
    })

  cat('<div id="bg">\n', paste(block, collapse = '\n'), '\n</div>')
}

# Bullet
bll = function(x) {
  paste0('* ', x, '\n')
  }

# Badge
badge = function(x) {
  paste0('[![Travis build status](https://travis-ci.org/caravagn/', x,'.svg?branch=master)](https://travis-ci.org/caravagn/', x, ')\n')
}

# entry
generator_entry = function(PKG)
{
  block = sapply(
    PKG,
    function(x){
      paste(
        paste0('\n-----\n\n **', x, ' - ', TPKG[[x]], '**'),
        isrcr(x), '\n\n',
        badge(x), '\n\n',
        bll(DPKG[[x]])
        )
    })

  paste('<div id="bg">\n', paste(block, collapse = '\n'), '\n</div>')
}

clipr::write_clip(generator_entry(PKG))

# install.packages('clipr')


