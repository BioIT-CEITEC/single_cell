##########################################
# wrapper for rule: seurat_object_creation
##########################################
import os
from snakemake.shell import shell
from os.path import dirname
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: seurat_object_creation \n##\n")
f.close()

command = "Rscript " + os.path.abspath(dirname(__file__) + "/seurat_obj.R ") + \
            dirname(snakemake.input.counts) + " " + dirname(snakemake.output.rds) + " "+ snakemake.params.sample  + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)