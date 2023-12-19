readextract4QC <- function(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/virusoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/Viral-Track/Virusite_annotation_file.txt",
                           readnameFilter = NULL){
  
  message("info::scanning files")
  library(dplyr)
  if (dir.exists(root_dir) == FALSE){
    stop("the root dir does not exist ")
  }
  sampledirs <-  fs::dir_ls(root_dir)
  detectstatus <- list()
  for ( i in 1:length(sampledirs)){
    chr_target <- readr::read_csv(glue::glue("{sampledirs[i]}/chr_target.txt"))
    chr_target <- chr_target$chr
    chr_name <- readr::read_csv(glue::glue("{sampledirs[i]}/chr_names.txt"))
    chr_name <- chr_name$chr
    if (!is.null(readnameFilter)){
      chr_target <- chr_target[grepl(readnameFilter,chr_target)]
      chr_name <- chr_name[grepl(readnameFilter,chr_name)]
      for (i in 1:length(chr_name)){
        chr_name[i] <- paste0(chr_name[i],".bam")
      }
      readr::write_csv(data.frame(chr = chr_target),
                       glue::glue("{sampledirs[i]}/chr_target.txt"),quote = "none",
                       col_names = FALSE)
      readr::write_csv(data.frame(chr = chr_name),
                       glue::glue("{sampledirs[i]}/chr_name.txt"),quote = "none",
                       col_names = FALSE)
    }
    output <- glue::glue("{sampledirs[i]}/{output_dir}")
    fs::dir_create(output)
  }
}

viraldetect <- readextract4QC(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/virusoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/Viral-Track/utils/Virusite_annotation_file.txt",
                           readnameFilter = "refseq")

