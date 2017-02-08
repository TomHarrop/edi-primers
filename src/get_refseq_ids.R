#!/usr/bin/env Rscript

library(data.table)

gene.ids <- fread("data/gene_list.txt")

# fix messy LOC IDs
gene.ids[, MsuId := gsub("Loc", "LOC", MsuId), by = MsuId]

refseq_list <- oryzr::LocToRefSeq(gene.ids[, unique(MsuId)])

# manually add records
refseq_list$notMatched[ensembl_gene_id == "OS02G0747600",
                       refseq_mrna := "XM_015770945"]
refseq_list$notMatched[ensembl_gene_id == "OS03G0186600",
                       refseq_mrna := "XM_015775124"]
refseq_list$notMatched[ensembl_gene_id == "OS03G0680200",
                       refseq_mrna := "XM_015774593"]
refseq_list$notMatched[ensembl_gene_id == "OS03G0770700",
                       refseq_mrna := "XM_015774163"]
refseq_list$notMatched[ensembl_gene_id == "OS04G0494950",
                       refseq_mrna := "XM_015779643"]
refseq_list$notMatched[ensembl_gene_id == "OS10G0536100",
                       refseq_mrna := "XM_015759313"]
refseq_list$notMatched[ensembl_gene_id == "OS04G0580700",
                       refseq_mrna := "XM_015778933"]
refseq_list$notMatched[ensembl_gene_id == "OS05G0118700",
                       refseq_mrna := "XM_015783597"]
refseq_list$notMatched[ensembl_gene_id == "OS05G0497200",
                       refseq_mrna := "XM_015782903"]
refseq_list$notMatched[ensembl_gene_id == "OS05G0567600",
                       refseq_mrna := "XM_015783987"]
refseq_list$notMatched[ensembl_gene_id == "OS06G0198500",
                       refseq_mrna := "XM_015786423"]
refseq_list$notMatched[ensembl_gene_id == "OS06G0691100",
                       refseq_mrna := "XM_015785536"]
refseq_list$notMatched[ensembl_gene_id == "OS08G0537900",
                       refseq_mrna := "XM_015794616"]
refseq_list$notMatched[ensembl_gene_id == "OS12G0617000",
                       refseq_mrna := "NM_001073818"]
refseq_list$notMatched[ensembl_gene_id == "OS02G0493050",
                       refseq_mrna := "XM_015768952"]
refseq_list$notMatched[tigrId == "LOC_Os09g20350",
                       refseq_mrna := "XM_015756652"]

refseq <- rbindlist(refseq_list, idcol = FALSE)
msu_refseq <- merge(gene.ids,
      refseq[, .(tigrId, refseq_mrna)],
      by.x = "MsuId",
      by.y = "tigrId",
      all.x = TRUE)

# write output for primer script
msu_refseq[refseq_mrna == "", refseq_mrna := NA]
msu_refseq[!is.na(refseq_mrna),
           write.table(unique(refseq_mrna),
                       file = "data/refseq.txt",
                       quote = FALSE,
                       row.names = FALSE,
                       col.names = FALSE)]
saveRDS(msu_refseq, 'output/msu_refseq.Rds')

