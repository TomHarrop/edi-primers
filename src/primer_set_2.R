library(data.table)

bof <- fread("data/primers_2017-02-08.csv")

# get an exhaustive list of RAP IDs
refseq_list <- oryzr::LocToRefSeq(bof[, unique(gene)])
refseq_results <- rbindlist(refseq_list, idcol = FALSE)
refseq_results[is.na(ensembl_gene_id) | ensembl_gene_id == "",
               ensembl_gene_id := NA]
tigr_to_rap_biomart <- refseq_results[,.(
    ensembl_gene_id = unlist(strsplit(ensembl_gene_id, "/", fixed = TRUE))),
    by = tigrId]
setnames(tigr_to_rap_biomart, 'tigrId', 'gene')
msu_to_rap_offline <- oryzr::LocToGeneName(bof[, unique(gene)])[, .(
    gene = MsuID, ensembl_gene_id = toupper(RapID))]

all_raps <- unique(merge(tigr_to_rap_biomart[!is.na(ensembl_gene_id)],
                         msu_to_rap_offline[!is.na(ensembl_gene_id)],
                         all = TRUE))

all_raps[, write.table(unique(ensembl_gene_id),
                       file = "test/rap_ids.csv",
                       sep = ",",
                       quote = FALSE,
                       na = "",
                       row.names = FALSE, col.names = FALSE)]


# now run rap_to_ncbi.py

# read the results back in
rap_to_gi <- fread("test/rap_to_gi.csv")
gi_to_nm <- fread("test/gi_to_nm.csv")
gi_to_xm <- fread("test/gi_to_xm.csv")

nm_merged <- merge(rap_to_gi, gi_to_nm, all = TRUE)
xm_merged <- merge(rap_to_gi, gi_to_xm, all = TRUE)

both_tx <- merge(nm_merged, xm_merged, by = c("gi", "rap"), all = TRUE)
both_tx[transcript.x == "" & transcript.y != "", transcript := transcript.y]
both_tx[transcript.y == "" & transcript.x != "", transcript := transcript.x]
both_tx[transcript.x != "" & transcript.y != ""]

both_tx[, c("transcript.x", "transcript.y") := NULL]

# merge to LOC table
all_genes_with_tx <- merge(all_raps, both_tx,
                           by.x = "ensembl_gene_id",
                           by.y = "rap", all = TRUE)
setkey(all_genes_with_tx, gene, transcript)
unique(all_genes_with_tx)

# now find missing transcripts, neaten table and run primer design

