######################################
# wrapper for rule: fastq_symlink
######################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_symlink \n##\n")
f.close()


# CREATE the symbolic link between the file in raw_fastq and singleCell_fastq folder
command = "mkdir -p " + os.path.dirname(snakemake.output.singleCell)
shell(command)

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

command = "cp " + snakemake.input.fastq +" "+snakemake.output.singleCell
shell(command)

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

