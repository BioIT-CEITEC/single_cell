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

#working_folder = os.getcwd() + "/results/" + snakemake.wildcards.sample + "/" + snakemake.wildcards.sample + "_Log.out"

command = "Rscript " + os.path.abspath(dirname(__file__) + "/seurat_obj.R ") + \
            dirname(snakemake.input.counts) + " " + dirname(snakemake.output.rds) + " "+ snakemake.wildcards.sample  + " >> " + log_filename + " 2>&1" #+ " " + dirname(snakemake.output.folder)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
