readextract4QC <- function(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/germoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/Viral-Track/Virusite_annotation_file.txt",
                           readnameFilter = NULL,
                           N_thread = 32){
  
  message("info::loading functions")
  ##Function to extract information from STAR Log output
  
  Extraction_Log_final = function(path_to_Log_file) {
    #Loading the final log file from STAR
    Log_file =read.delim(path_to_Log_file,header = F,sep="\t")
    Name_variable = as.character(Log_file$V1)
    Value_variable = as.character(Log_file$V2)
    
    
    #Extracting the informations about the mapping quality 
    Uniquely_mapped_percent = Value_variable[Name_variable=="                        Uniquely mapped reads % |"]
    Multiple_mapped_percent = Value_variable[Name_variable=="             % of reads mapped to multiple loci |"]
    Unmapped_too_short_percent = Value_variable[Name_variable=="                 % of reads unmapped: too short |"]
    Unmapped_mismatch_percent = Value_variable[Name_variable=="       % of reads unmapped: too many mismatches |"]
    Unmapped_other_percent = Value_variable[Name_variable=="                     % of reads unmapped: other |"]
    
    remove_percent = function(x) {
      l = nchar(x)
      x = substr(x,1,l-1)
      x=as.numeric(x)
      return(x)
    }  
    
    Uniquely_mapped_percent=remove_percent(Uniquely_mapped_percent)
    Multiple_mapped_percent=remove_percent(Multiple_mapped_percent)
    Unmapped_too_short_percent=remove_percent(Unmapped_too_short_percent)
    Unmapped_mismatch_percent = remove_percent(Unmapped_mismatch_percent)
    Unmapped_other_percent = remove_percent(Unmapped_other_percent)
    Total_unmapped = Unmapped_too_short_percent + Unmapped_mismatch_percent + Unmapped_other_percent
    
    Matrix_mapping = matrix(c(Uniquely_mapped_percent,Multiple_mapped_percent,Total_unmapped),nrow = 3)
    
    #####Looking at length of mappind, deletion and insertion
    
    Mean_mapped_length = as.numeric(Value_variable[Name_variable=="                          Average mapped length |"])
    Mean_deletion_length = as.numeric(Value_variable[Name_variable=="                        Deletion average length |"])
    Mean_insertion_length = as.numeric(Value_variable[Name_variable=="                       Insertion average length |"])
    #####Looking at rates of  of mismatch, insertion and deletion
    
    Mismatch_rate= remove_percent(Value_variable[Name_variable=="                      Mismatch rate per base, % |"])
    Deletion_rate = remove_percent(Value_variable[Name_variable=="                         Deletion rate per base |"])
    Insertion_rate = remove_percent(Value_variable[Name_variable=="                        Insertion rate per base |"])
    
    List_elements = list(Mapping_result= Matrix_mapping , 
                         Length_vector = c(Mean_mapped_length,Mean_insertion_length,Mean_deletion_length),
                         Rate_vector = c(Mismatch_rate,Insertion_rate,Deletion_rate) )
    
    return(List_elements)
  }
  
  
  ##Function to compute similarity between genome sequences
  
  alignement_score <- function(x) { # x is a vector of sequences in the Biostring format (DNAstring format)
    dist_matrix = matrix(0,nrow = length(x),ncol = length(x))
    comparison_sequence = c()
    for (i in 1:nrow(dist_matrix)){
      for (j in 1:ncol(dist_matrix)) {
        dist_matrix[i,j] = pairwiseAlignment(pattern = x[i], subject = x[j], type = "local-global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA") , gapOpening = 1, gapExtension = 0, scoreOnly=TRUE)
      }
    }
    dist_matrix = -dist_matrix
    rownames(dist_matrix) = names(x)
    colnames(dist_matrix) = names(x)
    return(dist_matrix)
  }
  
  string.to.colors = function (string, colors = NULL) 
  {
    if (is.factor(string)) {
      string = as.character(string)
    }
    if (!is.null(colors)) {
      if (length(colors) != length(unique(string))) {
        (break)("The number of colors must be equal to the number of unique elements.")
      }
      else {
        conv = cbind(unique(string), colors)
      }
    }
    else {
      conv = cbind(unique(string), rainbow(length(unique(string))))
    }
    unlist(lapply(string, FUN = function(x) {
      conv[which(conv[, 1] == x), 2]
    }))
  }
  
  #Loading of the libraries
  
  cat("Loading of the libraries.... ")
  
  suppressMessages(library(Biostrings))
  suppressMessages(library(ShortRead))
  suppressMessages(library(doParallel))
  suppressMessages(library(GenomicAlignments))
  cat("... done ! ")
  
  ##Registering the parallel environment
  cl =makeCluster(N_thread)
  registerDoParallel(cl)
  message("info::scanning files")
  library(dplyr)
  if (dir.exists(root_dir) == FALSE){
    stop("the root dir does not exist ")
  }
  Virus_database = read.delim(readannotation,header=T,sep="\t")
  sampledirs <-  fs::dir_ls(root_dir)
  detectstatus <- list()
  for ( i in 1:length(sampledirs)){
    dirid <- sampledirs[i]
    setwd(dirid)
    if (file.exists(glue::glue("{sampledirs[i]}/QCend.txt"))){
      next
    } else {
    bamfile <- list.files(sampledirs[i],"Aligned.sortedByCoord.out.bam",full.names = TRUE)
    bamfile <- bamfile[!grepl(".bam.bai",bamfile)]
    chr_target <- readr::read_csv(glue::glue("{sampledirs[i]}/chr_target.txt"))
    chr_target <- chr_target$chr
    message(glue::glue("info:: processing {i}/{length(sampledirs)} samples,working with {sampledirs[i]}"))
    if (!is.null(readnameFilter)){
      chr_target <- chr_target[!grepl("chr",chr_target)]
      chr_target <- chr_target[!grepl("NA/",chr_target)]
    }
    output <- glue::glue("{sampledirs[i]}/{output_dir}")
    outbam <- rep(NA,length(chr_target))
    for (i in 1:length(chr_target)){
      outbamfile <- glue::glue("{output}/{chr_target[i]}")
      outbam[i] <- outbamfile
    }
    
    message("calculating QC_result")
    
    print(outbam[1:6])
    QC_result = foreach(i=outbam,.combine = rbind,.packages = c("GenomicAlignments","ShortRead")) %dopar% {
      BAM_file= readGAlignments(i,param = ScanBamParam(what =scanBamWhat()))
      #Let's check the diversity of the reads
      Viral_reads = unique(BAM_file@elementMetadata$seq)
      Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob =T )
      Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]
      
      if ("numeric" %in% class(Viral_reads_contents)) {
        Viral_reads_contents = matrix(Viral_reads_contents_mean,ncol = 4)
      }
      
      Viral_reads_contents_mean =colMeans(Viral_reads_contents)
      Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)
      
      #... the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
      Covered_genome = coverage(BAM_file)[[i]]
      Covered_genome = as.numeric(Covered_genome)
      Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
      Covered_genome = rle(sign(Covered_genome))
      Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
      
      ##... the mean reads quality
      Reads_quality = as.character(BAM_file@elementMetadata$qual)
      Reads_quality = PhredQuality(Reads_quality)
      Reads_quality = as(Reads_quality,"IntegerList")
      Reads_quality = as.numeric(as.matrix(Reads_quality))
      Mean_read_quality = mean(Reads_quality)
      Sd_read_quality = sd(Reads_quality)
      
      ##... the number of mapped reads and unique mapped reads
      N_unique_mapped_reads = sum(BAM_file@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
      N_mapped_reads = length(BAM_file)
      Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
      
      
      ##... SDUST score to compute more efficiently the quality of the reads
      Mean_dust_score = NA
      Percent_high_quality_reads = NA
      if ("ShortRead"%in%installed.packages()){
        DUST_score = dustyScore( BAM_file@elementMetadata$seq)
        Mean_dust_score = mean(DUST_score)
        Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
      }
      
      ##... and lastly the pattern of the mapping 
      
      QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
                  Mean_read_quality,Sd_read_quality,
                  Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
                  Mean_dust_score,Percent_high_quality_reads)
      QC_temp
    }
    
    colnames(QC_result) = c("N_reads","N_unique_reads","Percent_uniquely_mapped",
                            "Mean_read_quality","Sd_read_quality",
                            c("A","C","G","T"),"Sequence_entropy","Spatial_distribution","Longest_contig",
                            "DUST_score","Percent_high_quality_reads")
    rownames(QC_result) = chr_target
    QC_result = as.data.frame(QC_result)
    QC_result = QC_result[QC_result$N_unique_reads>0,]
    message(paste("Exporting QC ",sampledirs[i],"...."))
    logfile <- list.files(getwd(),"Log.final.out")
    print(logfile)
    path_to_Log_file = paste(getwd(),"/",logfile,sep = "")
    Mapping_information = Extraction_Log_final(path_to_Log_file)
    Mean_mapping_length = Mapping_information$Length_vector[1]
    
    detected_virus = rownames(QC_result[QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])
    
    if (length(detected_virus)==0) {
      cat("No viral sequences detected detected \n")
    }
    
    if (length(detected_virus)>0) {
      cat(paste(length(detected_virus)," viral sequences detected detected \n",sep = ""))
    }
    
    cat(paste("Exporting QC table for ",dirid,"...."))
    
    #Removing low quality virus
    
    Filtered_QC=QC_result[detected_virus,]
    
    ##Exporting the tables of the QC analysis
    write.table(file = paste(dirid,"/","QC_unfiltered.txt",sep = ""),
                x = QC_result,quote = F,sep = "\t")
    write.table(file = paste(dirid,"/","QC_filtered.txt",sep = ""),
                x = Filtered_QC,quote = F,sep = "\t")
    cat(" done !")
    
    
    
    ##But also on the splicing events identified by STAR....
    
    Splice_table_path = list.files(dirid,
                                   pattern = "SJ.out.tab",full.names = T)
    Splice_table = read.table(Splice_table_path,header = F,sep = "\t")
    Splice_table = Splice_table[as.character(Splice_table$V1)%in%detected_virus,]
    colnames(Splice_table) = c("Virus","Start","End","Strand","Motif","Annotated","Uniquely mapped reads","Reads crossing the junction","Unknown")
    Splice_table$Annotated[Splice_table$Annotated==1]="GT/AG"
    Splice_table$Annotated[Splice_table$Annotated==2]="CT/AC"
    Splice_table$Annotated[Splice_table$Annotated==3]="GC/AG"
    Splice_table$Annotated[Splice_table$Annotated==4]="CT/GC"
    Splice_table$Annotated[Splice_table$Annotated==5]="AT/AC"
    Splice_table$Annotated[Splice_table$Annotated==6]="GT/AT"
    Splice_table$Annotated[Splice_table$Annotated==0]="Non-canonical"
    Splice_table$Annotated = factor(Splice_table$Annotated,levels = c("GT/AG","CT/AC","GC/AG","CT/GC","AT/AC","GT/AT","Non-canonical"))
    
    ###Additional info : % of reads mapped to viral vs host
    temp_chromosome_count_path <- paste0(dirid,"/Count_chromosomes.txt")
    Read_count_temp = read.table(temp_chromosome_count_path,header = F,row.names = 1)
    colnames(Read_count_temp) = c("Chromosome_length","Mapped_reads","Unknown")
    Read_count_temp = Read_count_temp[Read_count_temp$Mapped_reads!=0,]
    host_mapping_count = sum(Read_count_temp[grepl(pattern = "chr",rownames(Read_count_temp)),"Mapped_reads"])
    viral_mapping_count = sum(Read_count_temp[grepl(pattern = "NC",rownames(Read_count_temp)),"Mapped_reads"])
    total_mapping = viral_mapping_count + host_mapping_count
    Ratio_host_virus = matrix(data = c(host_mapping_count,viral_mapping_count)/total_mapping,ncol = 1)*100
    
    
    ##and also number of uniquely mapped reads and other reads for the filtered virus 
    
    Mapping_selected_virus = data.frame(Unique_mapping = (Filtered_QC$N_unique_reads),All_mapping = (Filtered_QC$N_reads),row.names = rownames(Filtered_QC))
    Mapping_selected_virus = Mapping_selected_virus[order(Mapping_selected_virus$Unique_mapping,decreasing = T),]
    
    
    ###Starting to plot the pdf QC
    if (length(detected_virus) > 0) {
      cat(paste("Creating QC plot for ",dirid,"...."))
      QClist <- list(
        Mapping_information,
        QC_result,
        Mapping_selected_virus
      )
      saveRDS(QClist,paste0(dirid,"/QC_report.RDS"))
    }
    ##Merging all viral sam files corresponding to identified viruses 
    cat("Merging Viral SAM files identified")
    
    list_BAM_files = paste(dirid,"/",output_dir,"/",sep="")
    selected_virus = list.files(list_BAM_files,full.names=F)
    selected_virus = base::strsplit(x = selected_virus,split = ".bam")
    selected_virus = unlist(lapply(selected_virus, function(x) {x[1]}))
    list_BAM_files = list.files(list_BAM_files,full.names=T)
    names(list_BAM_files) = selected_virus
    #filterqcvirus <- stringr::str_replace_all(rownames(Filtered_QC),"\\|","_")
    list_BAM_files = list_BAM_files[rownames(Filtered_QC)]
    list_BAM_files = paste("\'",list_BAM_files,"\'",sep = "")
    print(list_BAM_files)
    Merging_BAM_commad = paste(samtools_exe, "merge",paste(dirid,"Merged_viral_mapping.bam",sep = ""),list_BAM_files)
    if (length(detected_virus) > 0) {
      system(Merging_BAM_commad)
      dffinal <- data.frame(
        sample = dirid,
        status = "detect"
      )
      
    } else {
      dffinal <- data.frame(
        sample = dirid,
        status = "undetect"
      )
    }
    cat("Viral detection step done !")
    fs::file_create(glue::glue("{dirid}/QCend.txt"))
    detectstatus[[i]] <- dffinal
    }
  }
 return(detectstatus)
}

viraldetect <- readextract4QC(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/germoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/referG/all_contig_anno.txt",
                           readnameFilter = "refseq",
                           N_thread = 8)

