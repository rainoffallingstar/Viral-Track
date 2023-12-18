readextract4QC <- function(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/virusoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/Viral-Track/Virusite_annotation_file.txt",
                           readnameFilter = NULL){
  message("info::scanning files")
  if (dir.exists(root_dir) == FALSE){
    stop("the root dir does not exist ")
  }
  Virus_database = read.delim(readannotation,header=T,sep="\t")
  sampledirs <-  fs::dir_ls(root_dir)
  detectstatus <- list()
  for ( i in 1:length(sampledirs)){
    bamfile <- list.files(sampledirs[i],"Aligned.sortedByCoord.out.bam")
    chr_target <- readr::read_csv(glue::glue("{sampledirs[i]}/chr_target.txt"))
    chr_target <- chr_target$chr
    chr_name <- readr::read_csv(glue::glue("{sampledirs[i]}/chr_name.txt"))
    chr_name <- chr_name$chr
    if (!is.null(readnameFilter)){
      chr_target <- chr_target[grepl(readnameFilter,chr_target)]
      chr_name <- chr_name[grepl(readnameFilter,chr_name)]
    }
    output <- glue::glue("{sampledirs[i]}/{output_dir}")
    fs::dir_create(output)
    message(glue::glue("info:: start processing sample : {bamfile}, samtools extracting "))
    outbam <- rep(NA,length(chr_target))
    for (i in 1:length(chr_target)){
      outbamfile <- glue::glue("{output}/{chr_name[i]}.bam")
      #temp_export_bam_command = paste("source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh &&  conda activate py39 && conda run samtools view -b",temp_sorted_bam,i,">",paste(k,"Viral_BAM_files/",i,".bam ",sep = ""))
      temp_export_bam_command = glue::glue("{samtools_exe} view -b {bamfile} '{chr_target[i]}' > {outbamfile}")
      system(temp_export_bam_command)
      outbam[i] <- outbamfile
    }
    
    message("calculating QC_result")
    
    outbam <- outbam %>% na.omit()
    QC_result = foreach(i=outbam,.combine = rbind,.packages = c("GenomicAlignments","ShortRead")) %dopar% {
      BAM_file= readGAlignments(paste0(output,"/",i),param = ScanBamParam(what =scanBamWhat()))
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
    logfile <- list.files(sampledirs[i],"Log.final.out")
    path_to_Log_file = paste(sampledirs[i],"/",logfile,sep = "")
    Mapping_information = Extraction_Log_final(path_to_Log_file)
    Mean_mapping_length = Mapping_information$Length_vector[1]
    
    detected_virus = rownames(QC_result[QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])
    
    if (length(detected_virus)==0) {
      cat("No viral sequences detected detected \n")
    }
    
    if (length(detected_virus)>0) {
      cat(paste(length(detected_virus)," viral sequences detected detected \n",sep = ""))
    }
    
    cat(paste("Exporting QC table for ",name_target,"...."))
    
    #Removing low quality virus
    
    Filtered_QC=QC_result[detected_virus,]
    
    ##Exporting the tables of the QC analysis
    write.table(file = paste(sampledirs[i],"/","QC_unfiltered.txt",sep = ""),
                x = QC_result,quote = F,sep = "\t")
    write.table(file = paste(sampledirs[i],"/","QC_filtered.txt",sep = ""),
                x = Filtered_QC,quote = F,sep = "\t")
    cat(" done !")
    
    
    
    ##But also on the splicing events identified by STAR....
    
    Splice_table_path = list.files(sampledirs[i],
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
    temp_chromosome_count_path <- paste0(sampledirs[i],"/Count_chromosomes.txt")
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
    
    cat(paste("Creating QC plot for ",sampledirs[i],"...."))
    
    
    pdf(paste(sampledirs[i],"QC_report.pdf",sep = ""),height = 18,width = 12)
    par(las=1,mfrow=c(4,3),mar=c(6,6,6,4))
    Color_vector = c("lightskyblue1","orange","grey80")
    #Plotting the proportion of uniquely mapped reas, unmapped etc...
    barplot(Mapping_information$Mapping_result,ylim=c(0,100),xlim=c(0,5),ylab="Percentage of reads (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5)
    legend(x = 1.5,y=50,legend = c("Unmapped","Mapped to multiple loci","Uniquely mapped"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
    
    #Size of mapping, insertion and deletion
    barplot(Mapping_information$Length_vector,col="black",names.arg = c("Mapping length","Insertion length","Deletion length"),
            horiz = T,xlim=c(0,max(Mapping_information$Length_vector[1])*1.2),xlab="Nucleotide length",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
    
    #Rate of mismatch, deletion and insertion
    barplot(Mapping_information$Rate_vector,col="black",names.arg = c("Mismatch rate","Insertion rate","Deletion rate"),
            horiz = T,xlim=c(0,max(Mapping_information$Rate_vector[1])*1.2),xlab="Rate (%)",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
    
    #Ratio 
    Color_vector = c("darkred","grey")
    barplot(Ratio_host_virus,ylim=c(0,100),xlim=c(0,5),ylab="Mapping events (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5)
    legend(x = 1.5,y=50,legend = c("Viral mapping","Host mapping"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
    
    #First QC for the viral hits 
    
    if (length(detected_virus) == 0) {
      Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange")) ##Viral sequences that passed QC : green
    }
    
    if (length(detected_virus) > 0) {
      Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange","green")) ##Viral sequences that passed QC : green
    }
    
    
    plot(QC_result$N_unique_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",
         cex=1.5,xlab="N unique reads mapped",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="")
    abline(h=10,lwd=2,lty=2,col="grey")
    abline(v=Minimal_read_mapped,lwd=2,lty=2,col="grey")
    
    #Second QC for the viral hits 
    plot(QC_result$Sequence_entropy,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,
         cex=1.5,xlab="Sequence complexity",ylab="% Mapped genome",cex.lab=1.5,ylim=c(0,100),main="")
    abline(h=10,lwd=2,lty=2,col="grey")
    abline(v=1.2,lwd=2,lty=2,col="grey")
    
    #Third QC for the viral hits 
    
    plot(QC_result$Longest_contig,QC_result$DUST_score,pch=21,bg=Color_vector,
         cex=1.5,xlab="Longest contig (nt)",ylab="DUST score",cex.lab=1.4,main="")
    abline(v=3*Mean_mapping_length,lwd=2,lty=2,col="grey")
    
    #Number of reads for each filtered virus
    
    if (length(detected_virus) > 0) {
      
      barplot(Mapping_selected_virus$All_mapping[nrow(Mapping_selected_virus):1],
              col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
              xlim=c(0,max(Mapping_selected_virus$All_mapping)*1.2),xlab="Number of mapped reads",
              names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
      
      
      barplot(Mapping_selected_virus$Unique_mapping[nrow(Mapping_selected_virus):1],
              col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
              xlim=c(0,max(Mapping_selected_virus$Unique_mapping)*1.2),xlab="Number of uniquely mapped reads",
              names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
      
      
      #What are the identified viruses ?
      Virus_database_filtered = Virus_database[Virus_database$Name_sequence%in%detected_virus,]
      N_virus_identified = table(as.character(Virus_database_filtered$Virus_name))
      
      if (nrow(N_virus_identified)>0) {
        plot(NULL,xlim=c(0,10),ylim=c(0,length(detected_virus)),xaxt="n",yaxt="n",xlab="",ylab="")
        text(x=rep(3.5,length(detected_virus)),y=1:length(detected_virus),cex=0.8,
             Virus_database_filtered$Virus_name)
        text(x=rep(8,length(detected_virus)),y=1:length(detected_virus),cex=0.8,
             Virus_database_filtered$Name_sequence)
        
      }
      
    }
    
    dev.off()
    cat("QC plot done ! \n")
    
    
    ##Merging all viral sam files corresponding to identified viruses 
    cat("Merging Viral SAM files identified")
    
    list_BAM_files = paste(sampledirs[i],"/",output_dir,"/",sep="")
    selected_virus = list.files(list_BAM_files,full.names=F)
    selected_virus = base::strsplit(x = selected_virus,split = ".bam")
    selected_virus = unlist(lapply(selected_virus, function(x) {x[1]}))
    list_BAM_files = list.files(list_BAM_files,full.names=T)
    names(list_BAM_files) = selected_virus
    filterqcvirus <- stringr::str_replace_all(rownames(Filtered_QC),"\\|","_")
    list_BAM_files = list_BAM_files[rownames(filterqcvirus)]
    list_BAM_files = paste("\'",list_BAM_files,"\'",sep = "")
    
    Merging_BAM_commad = paste(samtools_exe, "merge",paste(sampledirs[i],"Merged_viral_mapping.bam",sep = ""),list_BAM_files)
    if (length(list_BAM_files) > 0) {
      system(Merging_BAM_commad)
      dffinal <- data.frame(
        sample = sampledirs[i],
        status = "detect"
      )
    } else {
      dffinal <- data.frame(
        sample = sampledirs[i],
        status = "undetect"
      )
    }
    cat("Viral detection step done !")
    detectstatus[[i]] <- dffinal
  }
 return(detectstatus)
}

viraldetect <- readextract4QC(samtools_exe = "source /public3/home/scg9946/miniconda3/etc/profile.d/conda.sh && conda run samtools",
                           root_dir = "/public3/home/scg9946/Virus-Track/virusoutput",
                           output_dir = "Viral_BAM_files",
                           readannotation = "/public3/home/scg9946/Viral-Track/utils/Virusite_annotation_file.txt",
                           readnameFilter = "refseq")

