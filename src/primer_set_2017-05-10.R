library(data.table)

# duplicate gene resolver
ResolveDuplicateGenes <- function(gene_id, transcipts_dt){
    
    my_dt <- copy(transcripts_dt)
    
    # remove NA transcripts for gene_id
    my_dt <- my_dt[!(gene == gene_id & is.na(transcript))]
    
    # find XM or NM transcripts
    my_transcripts <- my_dt[gene == gene_id, unique(transcript)]
    xm_tx <- grep("^XM", my_transcripts, value = TRUE)
    nm_tx <- grep("^NM", my_transcripts, value = TRUE)
    
    # prefer XM
    if (length(xm_tx > 0) & length(nm_tx > 0)) {
        return(my_dt[!transcript %in% nm_tx])
    }
    
    # otherwise we will have to search for primers for all transcripts :(
    my_dt
}

# read data
primer_list <- fread("data/primers_2017-05-10.csv",
                     header = FALSE,
                     col.names = "gene")
outdir <- paste("output", Sys.Date(), sep = "/")
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

# get an exhaustive list of RAP IDs
refseq_list <- oryzr::LocToRefSeq(primer_list[, unique(gene)])
refseq_results <- rbindlist(refseq_list, idcol = FALSE)
refseq_results[is.na(ensembl_gene_id) | ensembl_gene_id == "",
               ensembl_gene_id := NA]
refseq_results[is.na(refseq_mrna) | refseq_mrna == "",
               refseq_mrna := NA]
tigr_to_rap_biomart <- refseq_results[,.(
    ensembl_gene_id = unlist(strsplit(ensembl_gene_id, "/", fixed = TRUE))),
    by = tigrId]
setnames(tigr_to_rap_biomart, 'tigrId', 'gene')
msu_to_rap_offline <- oryzr::LocToGeneName(primer_list[, unique(gene)])[, .(
    gene = MsuID, ensembl_gene_id = toupper(RapID))]

all_raps <- unique(merge(tigr_to_rap_biomart[!is.na(ensembl_gene_id)],
                         msu_to_rap_offline[!is.na(ensembl_gene_id)],
                         all = TRUE))

all_raps[, write.table(unique(ensembl_gene_id),
                       file = paste(outdir, "rap_ids.csv", sep = "/"),
                       sep = ",",
                       quote = FALSE,
                       na = "",
                       row.names = FALSE, col.names = FALSE)]

# now run rap_to_transcript.py
# should make this callable from R with a tmp file

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

# find missing transcripts from original oryzr search
all_genes_with_tx[is.na(transcript),
                  transcript := refseq_results[tigrId == gene, refseq_mrna],
                  by = gene]

# remove duplicate rows
unique_transcripts <- all_genes_with_tx[!duplicated(
    all_genes_with_tx, by = c("gene", "transcript"))]


# write output
unique_transcripts[!is.na(transcript),
                   write.table(unique(transcript),
                               file = paste(outdir, "refseq.txt",
                                            sep = "/"),
                               sep = ",",
                               quote = FALSE,
                               na = "",
                               row.names = FALSE, col.names = FALSE)]

# combine results
primer_results <- fread("output/2017-02-10/primers.csv")
gene_primers <- merge(all_genes_with_tx, primer_results, all = TRUE, by.x = "transcript", by.y = "ref_seq")

gene_primers[, symbol := oryzr::LocToGeneName(gene)[, symbols], by = gene]
setcolorder(gene_primers,
            c("gene", "symbol", "transcript", "gi", "status", "F", "TM_F", "R", 
              "TM_R", "product_size", "intron_size", "ensembl_gene_id"))
setorder(gene_primers, symbol, gene, na.last = TRUE)
setnames(gene_primers, "transcript", "primer_search")
setnames(gene_primers, "ensembl_gene_id", "rap_id")


write.csv(
    gene_primers[!duplicated(gene_primers, by = c("gene", "symbol", "gi"))],
    file = "output/2017-02-10/annotated_primers.csv",
    quote = FALSE,
    na = "",
    row.names = FALSE)


