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

command = "Rscript -e \"rmarkdown::render('{0}', output_file = '../../{1}', params=list(args1='{2}', args2='{3}'))\" >> {4} 2>&1".format(os.path.abspath(os.path.dirname(__file__) + "/qual_report.Rmd"),snakemake.output.html, snakemake.input.csv, snakemake.params.sample, log_filename)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)