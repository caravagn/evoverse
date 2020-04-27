# VAF Adjustment for CN and purity as 1/2*CCF
vaf_adjustCN = function(v, m, M, p, mut.allele =1)
{
  CN = m+M

  0.5 * v * ((CN-2) * p + 2) / (mut.allele * p)
}
