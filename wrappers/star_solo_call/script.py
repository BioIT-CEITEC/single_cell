######################################
# wrapper for rule: STARsolo_mapping
######################################
import os
from snakemake.shell import shell
from os.path import dirname
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: STARsolo_mapping \n##\n")
f.close()

if snakemake.params.assay == "tenX_v3" or snakemake.params.assay == "tenX_multiome":
    flags = " --soloCBwhitelist " + os.path.join(snakemake.params.tool_folder, "cellranger/barcodes", "3M-february-2018.txt") + " --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4"
if snakemake.params.assay == "tenX_v2":
    flags = " --soloCBwhitelist " + os.path.join(snakemake.params.tool_folder, "cellranger/barcodes", "737K-august-2016.txt") + " --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4"
if snakemake.params.assay == "tenX_5p":
    flags = " --soloCBwhitelist " + os.path.join(snakemake.params.tool_folder, "cellranger/barcodes", "737K-arc-v1.txt") + " --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB"
if snakemake.params.assay == "tenX_5p_pe":
    flags = " --soloCBwhitelist " + os.path.join(snakemake.params.tool_folder, "cellranger/barcodes", "737K-august-2016.txt") + " --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloBarcodeMate 1 --clip5pNbases 39 0 --soloUMIdedup 1MM_CR"
if snakemake.params.assay == "DropSeq" or snakemake.params.assay == "SeqWell":
    flags = " --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8"
if snakemake.params.assay == "ShareSeq":
    flags = " --soloCBlen 24 --soloUMIstart 25 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4"

if snakemake.params.cellfilt == "EmptyDrops_CR":
    flags = flags + " --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 "
else:
    flags = flags + " --soloCellFilter CellRanger2.2 3000 0.99 10 "


command = "mkdir -p " + dirname(snakemake.params.prefix) + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "STAR --genomeDir " + dirname(snakemake.input.index) + \
          " --readFilesCommand zcat " + \
          " --runThreadN " + str(snakemake.threads) + \
          " --readFilesIn " + str(snakemake.input.c2) + " " + str(snakemake.input.c1) + \
          " --soloType CB_UMI_Simple " + \
          " --outFileNamePrefix " + snakemake.params.prefix +  \
          " --limitGenomeGenerateRAM " + str(snakemake.resources.mem * 1000000000) + \
          " --soloCBstart 1"  + flags + \
          " --soloFeatures Gene" + \
          " --outSAMtype BAM Unsorted" + \
          " --outFilterType BySJout" + \
          " --alignIntronMax 100000" + \
          " --quantMode GeneCounts" + \
          " --outFilterScoreMin 30" + \
          " --outSAMattributes NH HI nM AS GX GN sS sQ sM " + \
          " --outTmpDir " + snakemake.params.tmpd + "/" + snakemake.params.sample + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools sort -o " + str(snakemake.output.bam) + " " + str(snakemake.params.prefix) + "Aligned.out.bam " 
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)