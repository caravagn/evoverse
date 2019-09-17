#' Print for mvMOBSTER objects.
#'
#' @param x A mvMOBSTER \code{mbst_data} object.
#' @param ...
#'
#' @return Nothing
#'
#' @export
#'
#' @examples
#' TODO
print.mbst_data = function(x, ...)
{
  pio::pioHdr(
    paste("mvMOBSTER dataset"),
    toPrint =
      c(
        `  Dataset` = x$description,
        `  Samples` = paste(x$samples, collapse = ', '),
        `Mutations` = paste('N =', N(x))
      )
  )

  prt = x$data %>%
    group_by(sample, variable) %>%
    filter(value > 0) %>%
    summarize(
      mean = round(mean(value), 2),
      median = round(median(value), 2),
      min = round(min(value), 2),
      max = round(max(value), 2)
    ) %>%
    arrange(variable)

  prt = prt %>% mutate(combo =
                         paste0(
                           sprintf('%3s', min),
                           ' | ',
                           sprintf('%3s', median),
                           ' | ',
                           sprintf('%3s', max)
                         )) %>%
    select(-mean,-median,-min,-max) %>%
    ungroup()

  df = lapply(x$samples,
              function(s) {
                sp = prt %>%
                  filter(sample == s) %>% ungroup()

                sp %>% spread(variable, combo)
              })

  df = Reduce(bind_rows, df)
  # colnames(df)[2:4] = c(
  #   'DP min | median | max',
  #   'NV min | median | max',
  #   'VAF min | median | max')

  pio::pioTit('Coverage (DP), reads with mutant (NV) and allele frequency (VAF) -- min | median | max')
  print(data.frame(df, check.names = FALSE))



  # prt %>% ungroup() %>%
  #   filter(variable != 'VAF')
  #   # select(-mean, -min, -max) %>%
  #   spread(variable, c(median, min, max))

  # prt = prt %>%
  #   nest() %>%
  #   unlist(recursive = FALSE)
  #
  # prt
  #

  # cat(
  #   sprintf(
  #     ''
  #     )
  # )
  #
  # for(i in 1:nrow(prt))
  # {
  #   cat(prt$sample[i], "")
  #
  # }

  # prt = function(x,
  #                cols,
  #                v)
  # {
  #   d = x$data[, cols]
  #   f = apply(d, 2, summary)
  #   f = apply(f, 2, round, digit = 2)
  #
  #   colnames(f) = x$samples
  #
  #   f = f[c(1,3,4,6), , drop = FALSE]
  #
  #   frmt = function(f, v)
  #   {
  #     if(f < v[1]) cat(crayon::red(sprintf('%8s', f)))
  #     if(f >= v[1] & f <= v[2]) cat(crayon::yellow(sprintf('%8s', f)))
  #     if(f > v[2]) cat(crayon::green(sprintf('%8s', f)))
  #
  #   }
  #
  #   for(i in 1:nrow(f)) {
  #     cat(sprintf('  %10s', rownames(f)[i]))
  #     for(j in 1:ncol(f)) {
  #       frmt(f[i, j], v)
  #     }
  #
  #     cat('\n')
  #   }
  #
  # }
  #

  # print(prt)

  if (!is.null(x$annotation))
  {
    lbl = head(unique(x$annotation$variable))

    # pio::pioTit('ANNOTATION')
    # print((x$annotation))
    pio::pioStr('Annotations', prefix = '\n\n',
                paste0(c(lbl, '...'), collapse = ', '))
  }

  cat('\n')
  pio::pioTit('LOGGED OPERATIONS')

  print(x$operationLog)

}
