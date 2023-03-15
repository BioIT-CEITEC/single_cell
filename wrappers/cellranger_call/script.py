######################################
# wrapper for rule: cellranger_call
######################################
import os
import uuid
from snakemake.shell import shell

shell.executable("/bin/bash")

pwd = os.getcwd()
log_filename = os.path.join(pwd,str(snakemake.log))

SAMPLE = [x for x in snakemake.params.sample]
LIBS = list(set([x.rsplit("_", 1)[0] for x in SAMPLE]))

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cellranger_call \n##\n")
f.close()


# CREATE the csv file if not existing and add the header
libfile = os.path.join(pwd,snakemake.params.libraries)
f = open(log_filename, 'at')
f.write("## INFO: creating " + libfile + " file\n")
f.close()
lf = open(libfile, "w")
lf.write("fastqs,sample,library_type\n")
lf.close()

# ADD info to csv
# CHECK if "GE" is in the name of the library
for x in LIBS:

    f = open(log_filename, 'at')
    f.write("## INFO: creating " + x + " folder\n")
    f.close()

    SAMPLE_LIB_DIR = os.path.join(pwd,"singleCell_fastq", x)
    line_to_write = SAMPLE_LIB_DIR + "," + x + "," + snakemake.params.library_types_dict[x] + "\n"
    lf = open(libfile, "at")
    lf.write(line_to_write)
    lf.close()

f = open(log_filename, 'at')
f.write("## INFO: filling " + libfile + " file\n")
f.close()

# Extract cellranger tar archive into tmp
command = "tar -zxvf "+snakemake.input.tar+" -C "+snakemake.params.wdir+" --skip-old-files >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# CALL cellranger
if snakemake.params.sc_hashtags != "no":
    feature_ref_parameter = "--feature-ref=" + snakemake.params.feature_ref_dir + "/feature_ref_" + snakemake.params.sc_hashtags + ".csv"
else:
    feature_ref_parameter = ""

workdir = snakemake.params.outdir+"_"+uuid.uuid4().hex
command = "cd " + snakemake.params.wdir + \
          " ; export PATH="+snakemake.params.bin_path+":$PATH && $(which time) cellranger count " + \
          " --id=" + workdir + \
          " --libraries=" + libfile + \
          " " + feature_ref_parameter + \
          " --transcriptome=" + snakemake.params.transcriptome + \
          " --localcores=" + str(snakemake.threads) + \
          " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "rsync -rtv " + os.path.join(snakemake.params.wdir,workdir,".") + " " + os.path.join(pwd,snakemake.params.outdir)+ " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
