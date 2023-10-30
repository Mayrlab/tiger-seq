library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(plyranges)
library(tidyverse)
library(magrittr)

gr <- read_gff("data/gff/MANE.GRCh38.v1.2.refseq_genomic.gff.gz")

txdb <- makeTxDbFromGRanges(gr)

named_nucleotide_fraction <- function (x, rowname_key="gene") {
    alphabetFrequency(x)[,c("A","C","G","T")] %>%
        `colnames<-`(str_c("frac_", c("A","C","G","T"))) %>%
        { . / rowSums(.) } %>%
        `rownames<-`(names(x)) %>%
        as_tibble(rownames=rowname_key)
}

gr_mrna <- gr %>% filter(type == "mRNA") %>%
    group_by(gene) %>% filter(rank(-width, ties.method="first") == 1) %>%
    ungroup() %>% `names<-`(.$gene)

tx2gene <- setNames(names(gr_mrna), gr_mrna$Name)

df_fracs_tx_unit <- gr_mrna %>%
    getSeq(x=BSgenome.Hsapiens.UCSC.hg38) %>%
    named_nucleotide_fraction

df_fracs_mrna <- txdb %>%
    extractTranscriptSeqs(x=BSgenome.Hsapiens.UCSC.hg38, use.names=TRUE) %>%
    named_nucleotide_fraction(rowname_key="tx_name") %>%
    filter(tx_name %in% names(tx2gene)) %>%
    mutate(gene=tx2gene[tx_name], tx_name=NULL)

df_fracs <- left_join(df_fracs_tx_unit, df_fracs_mrna, by="gene",
                      suffix=c("_tx_unit", "_mRNA"))

df_genes <- gr %>%
    filter(type == "mRNA") %>%
    as_tibble() %>%
    group_by(gene) %>%
    summarize(n_txs=dplyr::n(),
              tx_unit_length=max(width)) %>%
    left_join(df_fracs, by="gene") %>%
    arrange(-n_txs)

df_in <- readxl::read_xlsx("data/CM table 18388 anno.xlsx")


gr %>%
    filter(type == "gene") %>%
    mutate(gene_id=str_replace(sapply(Dbxref, str_c, collapse=";"), "GeneID:([0-9]+);.*", "\\1"),
           hgnc_id=str_replace(sapply(Dbxref, str_c, collapse=";"), ".*HGNC:([0-9]+)(;.*)?", "\\1")) %>%
    filter(str_detect(gene_id, "44996") | str_detect(hgnc_id, "44996"))

gr_genes <- gr %>%
    filter(type == "gene") %>%
    mutate(gene_id=str_replace(sapply(Dbxref, str_c, collapse=";"), "GeneID:([0-9]+);.*", "\\1"),
           hgnc_id=str_replace(sapply(Dbxref, str_c, collapse=";"), ".*HGNC:([0-9]+)(;.*)?", "\\1"))


idx_gene <- df_in$`Gene...4` %in% gr_genes$gene

df_in[!idx_gene,]

df_exons_all <- gr %>%
    filter(type == "exon") %>%
    as_tibble() %>%
    mutate(Parent=unlist(Parent)) %>%
    arrange(ID) %>%
    group_by(gene, Parent) %>%
    summarize(mRNA_length=sum(width),
              n_exons_all=dplyr::n(),
              mean_exon_length_all=mean(width),
              first_exon_length=first(width),
              last_exon_length=last(width),
              .groups='drop')

df_exons_cds <- gr %>%
    filter(type == "CDS") %>%
    as_tibble() %>%
    mutate(Parent=unlist(Parent)) %>%
    arrange(ID) %>%
    group_by(gene, Parent) %>%
    summarize(cds_length=sum(width),
              n_exons_cds=dplyr::n(),
              mean_cds_exon_length_all=mean(width),
              first_cds_exon_length=first(width),
              last_cds_exon_length=last(width),
              .groups='drop') %>%
    mutate(mean_cds_exon_length_internal=ifelse(n_exons_cds < 3, NA,
                                                (cds_length - first_cds_exon_length - last_cds_exon_length)/(n_exons_cds - 2)))

df_exons_combined <- full_join(df_exons_all, df_exons_cds, by=c("gene", "Parent")) %>%
    group_by(gene) %>%
    slice_min(order_by=tibble(-mRNA_length, -n_exons_all, Parent)) %>%
    ungroup() %>%
    dplyr::rename(transcript_id_mane=Parent) %>%
    left_join(df_genes, by="gene")

genes_missing <- df_in$Gene...4[!(df_in$Gene...4 %in% df_exons_combined$gene)] %>% unique()

map_gene_syns <- gr_genes %>%
    as_tibble %>% select(gene, gene_synonym) %>%
    deframe()

res <- lapply(genes_missing, function (gene) {
    lapply(map_gene_syns, function (x) gene %in% x) %>% unlist %>% { names(.)[.] }
}) %>% `names<-`(genes_missing)


map2gene <- res %>%
    { .[lengths(.) > 0] } %>%
    { sapply(., first) }

df_out <- df_in %>%
    mutate(gene_key=ifelse(`Gene...4` %in% names(map2gene), map2gene[`Gene...4`], `Gene...4`),
           used_gene_synonym=gene_key != `Gene...4`) %>%
    left_join(df_exons_combined, by=c("gene_key"="gene"))

writexl::write_xlsx(df_out, path="out/gene_annots_mane.xlsx")

