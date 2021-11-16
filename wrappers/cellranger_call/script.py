######################################
# wrapper for rule: cellranger_call
######################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


SAMPLE = [x for x in snakemake.params.sample]
LIBS = list(set([x.rsplit("_", 1)[0] for x in SAMPLE]))

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cellranger_call \n##\n")
f.close()


# CREATE the csv file if not existing and add the header
lf = open(snakemake.params.libraries, "w")
lf.write("fastqs,sample,library_type\n")
lf.close()
f = open(log_filename, 'at')
f.write("## COMMAND: create" + snakemake.params.libraries + " file\n")
f.close()

# ADD info to csv
# CHECK if "GE" is in the name of the library
for x in LIBS:

    f = open(log_filename, 'at')
    f.write("## COMMAND: create " + x + " folder\n")
    f.close()

    SAMPLE_LIB_DIR = os.path.join(snakemake.params.wdir, "singleCell_fastq", x)
    line_to_write = SAMPLE_LIB_DIR + "," + x + "," + snakemake.params.library_types_dict[x] + "\n"
    lf = open(snakemake.params.libraries, "at")
    lf.write(line_to_write)
    lf.close()

f = open(log_filename, 'at')
f.write("## COMMAND: filling" + snakemake.params.libraries + " file\n")
f.close()

# CALL cellranger
if snakemake.params.sc_hashtags != "no":
    feature_ref_parameter = "--feature-ref=" + snakemake.params.feature_ref_dir + "/feature_ref_" + snakemake.params.sc_hashtags + ".csv"
else:
    feature_ref_parameter = ""

command = "cd " + snakemake.params.wdir + " ; rm -Rf " + snakemake.params.outdir + " ; " + snakemake.params.binary + " count " + \
          " --id=" + snakemake.params.outdir + \
          " --libraries=" + os.path.basename(snakemake.params.libraries) + \
          " " + feature_ref_parameter + \
          " --transcriptome=" + snakemake.params.transcriptome + \
          " --localcores " + str(snakemake.threads) + \
          " >> " + log_filename.replace(snakemake.params.wdir + "/", "") + " 2>&1 ; cd .."

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# command = "cp " + snakemake.params.HTML + " " + snakemake.output.REPORT
# f = open(log_filename, 'at')
# f.write("## COMMAND: " + command + "\n")
# f.close()
# shell(command)
