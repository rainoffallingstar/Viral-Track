# in virus mod
echo "in the virus mod"
echo "start scanning"
Rscript Viral_Track_scanning.R /public3/home/scg9946/Viral-Track/virus_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt
echo "Track_transcript_assembly"
Rscript Viral_Track_transcript_assembly.R /public3/home/scg9946/Viral-Track/virus_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt
echo "cell_demultiplexing"
Rscript Viral_Track_cell_demultiplexing.R /public3/home/scg9946/Viral-Track/virus_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt

# in germ mod 
echo "in the germ mod"
echo "start scanning"
Rscript Viral_Track_scanning.R /public3/home/scg9946/Viral-Track/germ_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt
echo "Track_transcript_assembly"
Rscript Viral_Track_transcript_assembly.R /public3/home/scg9946/Viral-Track/germ_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt
echo "cell_demultiplexing"
Rscript Viral_Track_cell_demultiplexing.R /public3/home/scg9946/Viral-Track/germ_Parameters.txt /public3/home/scg9946/Viral-Track/Files_to_process.txt