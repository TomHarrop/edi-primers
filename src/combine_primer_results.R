#!/usr/bin/env Rscript

library(data.table)

primers <- fread("output/primers.csv")
msu.refseq <- readRDS("output/msu_refseq.Rds")

primer.results <- merge(msu.refseq,
      primers,
      by.x = "refseq_mrna",
      by.y = "ref_seq",
      all = TRUE)

primer.results[, symbol := oryzr::LocToGeneName(MsuId)[, symbols], by = MsuId]
setcolorder(primer.results,
         c("MsuId", "symbol", "refseq_mrna", "status", "F", "TM_F", "R", 
           "TM_R", "product_size", "intron_size"))
setorder(primer.results, symbol, MsuId, na.last = TRUE)
setnames(primer.results, "refseq_mrna", "primer_search")

write.csv(primer.results,
          file = "output/annotated_primers.csv",
          quote = FALSE,
          na = "",
          row.names = FALSE)