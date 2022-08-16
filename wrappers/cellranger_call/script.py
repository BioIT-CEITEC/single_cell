######################################
# wrapper for rule: cellranger_call
######################################
import os, zipfile

from snakemake.shell import shell

shell.executable("/bin/bash")
log_filename = str(snakemake.log)


# unzip cellranger

def extract_file(zip_file, info, extract_dir):
    zip_file.extract(info.filename, path=extract_dir)
    out_path = os.path.join(extract_dir, info.filename)

    perm = info.external_attr >> 16
    os.chmod(out_path, perm)


with zipfile.ZipFile(snakemake.input.cellranger_dir_zip, 'r') as cellranger_zip:
    for info in cellranger_zip.infolist():
        extract_file(cellranger_zip, info, snakemake.params.wdir)

SAMPLE = [x for x in snakemake.params.sample]
LIBS = list(set([x.rsplit("_", 1)[0] for x in SAMPLE]))

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cellranger_call \n##\n")
f.close()

###### FASTQ SYMLINK preprocessing

# CREATE the symbolic link between the file in raw_fastq and singleCell_fastq folder
for singlecell_fastq in set(map(os.path.dirname, snakemake.output.c1 + snakemake.output.c2)):
    command = "mkdir -p " + singlecell_fastq
    shell(command)
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()


for i in range(len(snakemake.input.fastq1)):
    command = "cp " + snakemake.input.fastq1[i] + " " + snakemake.output.c1[i]
    shell(command)

for i in range(len(snakemake.input.fastq2)):
    command = "cp " + snakemake.input.fastq2[i] + " " + snakemake.output.c2[i]
    shell(command)

# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()

# CREATE the csv file if not existing and add the header
lf = open(snakemake.input.libraries, "w")
lf.write("fastqs,sample,library_type\n")
lf.close()
f = open(log_filename, 'at')
f.write("## COMMAND: create" + snakemake.input.libraries + " file\n")
f.close()

# ADD info to csv
# CHECK if "GE" is in the name of the library
for x in LIBS:
    f = open(log_filename, 'at')
    f.write("## COMMAND: create " + x + " folder\n")
    f.close()

    SAMPLE_LIB_DIR = os.path.join(snakemake.params.wdir, "singleCell_fastq", x)
    line_to_write = "/tmp/" + SAMPLE_LIB_DIR + "," + x + "," + snakemake.params.library_types_dict[x] + "\n"
    lf = open(snakemake.input.libraries, "at")
    lf.write(line_to_write)
    lf.close()

f = open(log_filename, 'at')
f.write("## COMMAND: filling" + snakemake.input.libraries + " file\n")
f.close()

# CALL cellranger
if snakemake.params.sc_hashtags != "no":
    feature_ref_parameter = "--feature-ref=" + "/tmp/" + snakemake.input.feature_ref_path
else:
    feature_ref_parameter = ""

cmd = "chmod +x " + snakemake.params.wdir + "/cellranger-5.0.1/bin/cellranger"
shell(cmd)

command = "cd " + snakemake.params.wdir + "; rm -Rf " + snakemake.params.outdir + " ; " + "cellranger-5.0.1/bin/cellranger" + " count " + \
          " --id=" + snakemake.params.outdir + \
          " --libraries=" + os.path.basename(snakemake.input.libraries) + \
          " " + feature_ref_parameter + \
          " --transcriptome=" + "/tmp/" + snakemake.params.transcriptome + \
          " --localcores " + str(snakemake.threads)
          # " >> " + log_filename.replace(snakemake.params.wdir + "/", "") + " 2>&1 ; cd .."

# command =  "cd " + snakemake.params.wdir + " ; rm -Rf " + snakemake.params.outdir + " ; " + snakemake.input.binary + " count " + \
#           " --id=" + snakemake.params.outdir + \
#           " --libraries=" + os.path.basename(snakemake.input.libraries) + \
#           " " + feature_ref_parameter + \
#           " --transcriptome=" + snakemake.input.transcriptome + \
#           " --localcores " + str(snakemake.threads) + \
#           " >> " + log_filename.replace(snakemake.params.wdir + "/", "") + " 2>&1 ; cd .."

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# command = "cp " + snakemake.params.HTML + " " + snakemake.output.REPORT
# f = open(log_filename, 'at')
# f.write("## COMMAND: " + command + "\n")
# f.close()
# shell(command)
