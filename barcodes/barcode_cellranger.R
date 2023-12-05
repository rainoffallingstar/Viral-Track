# processing barcodes from cellranger

library(readr)
barcodes_tsv <- read_csv("~/barcodes_tsv.txt", 
                         col_names = FALSE)
colnames(barcodes_tsv) <- "path_dir"

barcodes_tsv <- barcodes_tsv$path_dir[grepl("filtered_feature_bc_matrix",barcodes_tsv$path_dir)]
barcodes_tsv <- barcodes_tsv[grepl("./multi/pipeline/01_cellranger_multi/",barcodes_tsv)]
for (i in 1:length(barcodes_tsv)){
  barcodes_tsvtemp <- stringr::str_remove(barcodes_tsv[i],"\\.")
  barcodes_tsv[i] <- paste0("/public3/home/scg5695/project/15_scCRC",barcodes_tsvtemp)
}
library(dplyr)
barcode_df <- data.frame(
  path_dir = barcodes_tsv,
  ids = c(1:length(barcodes_tsv))
) %>% 
  dplyr::mutate(
    newpath = paste0("/public3/home/scg9946/barcodes/",ids,".tsv.gz")
  )
saveRDS(barcode_df,"barcode_df.RDS")
for (i in 1:nrow(barcode_df)){
  fs::dir_create("/public3/home/scg9946/barcodes/")
  fs::file_copy(barcode_df$path_dir[i],barcode_df$newpath[i])
}

barcode_list <- list()

for (i in 1:nrow(barcode_df)){
  barcode_list[[i]] <- read.table(glue::glue("barcodes/{i}.tsv.gz"))
}

barcodes <- data.table::rbindlist(barcode_list) %>% 
  as.data.frame() %>% 
  dplyr::distinct(V1)

barcodes <- purrr::map(barcodes$V1,function(x) stringr::str_remove(x,"-1")) %>% 
  unlist()
write_tsv(data.frame(barcodes = barcodes),"barcodes.txt",col_names = FALSE)

