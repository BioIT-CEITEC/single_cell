######################################
# wrapper for rule: cellranger_call
######################################
import os

from snakemake.shell import shell

# Computing environment preparation

shell.executable("/bin/bash")
log_filename = str(snakemake.log)

if snakemake.params.computing_type == "kubernetes":
    script_wdir = os.getcwd()

    # unzip cellranger
    command = "tar -xpf " + snakemake.input.cellranger_dir_tar + " -C " + snakemake.params.task_wdir
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    cellranger_binary = "cellranger-5.0.1/bin/cellranger"
else:
    script_wdir = ""
    cellranger_binary = snakemake.input.cellranger_dir_tar.replace(".tar.gz", "/bin/cellranger")

###### FASTQ preprocessing
# cellranger requires libraries to conform to a specific naming scheme

# create target directories
for singlecell_fastq in set(map(os.path.dirname, snakemake.params.c1 + snakemake.params.c2)):
    command = "mkdir -p " + singlecell_fastq
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)


# rename fastq
def prepare_fastq(raw_fastqs, output_fastqs):
    assert len(raw_fastqs) == len(output_fastqs)
    for i in range(len(raw_fastqs)):
        if snakemake.params.computing_type == "kubernetes":
            command = "mv " + raw_fastqs[i] + " " + output_fastqs[i]
        else:
            command = "ln -s " + os.path.join(snakemake.params.task_wdir, raw_fastqs[i]) + " " + os.path.join(snakemake.params.task_wdir, output_fastqs[i])
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)


prepare_fastq(snakemake.input.fastq1, snakemake.params.c1)
prepare_fastq(snakemake.input.fastq2, snakemake.params.c2)


# CREATE the csv file and add the header
lf = open(snakemake.params.libraries, "w")
lf.write("fastqs,sample,library_type\n")
lf.close()
f = open(log_filename, 'at')
f.write("## COMMAND: create " + snakemake.params.libraries + " file\n")
f.close()

# ADD info to csv
# CHECK if "GE" is in the name of the library
for x in snakemake.params.libs:
    f = open(log_filename, 'at')
    f.write("## COMMAND: create " + x + " folder\n")
    f.close()

    SAMPLE_LIB_DIR = os.path.join(script_wdir, snakemake.params.task_wdir, "singleCell_fastq", x)
    line_to_write = SAMPLE_LIB_DIR + "," + x + "," + snakemake.params.library_types_dict[x] + "\n"
    lf = open(snakemake.params.libraries, "at")
    lf.write(line_to_write)
    lf.close()

f = open(log_filename, 'at')
f.write("## COMMAND: filling  " + snakemake.params.libraries + " file\n")
f.close()

if snakemake.params.sc_hashtags != "no":
    feature_ref_parameter = "--feature-ref=" + os.path.join(script_wdir, snakemake.input.feature_ref_path)
else:
    feature_ref_parameter = ""

# CALL cellranger

command = "cd " + snakemake.params.task_wdir + "; rm -Rf " + snakemake.params.outdir + " ; " + cellranger_binary + " count " + \
          " --id=" + snakemake.params.outdir + \
          " --libraries=" + os.path.basename(snakemake.params.libraries) + \
          " " + feature_ref_parameter + \
          " --transcriptome=" + os.path.join(script_wdir, snakemake.params.transcriptome_files_path) + \
          " --localcores " + str(snakemake.threads) + \
          " >> " + log_filename.replace(snakemake.params.task_wdir + "/", "") + " 2>&1 ; cd .."


f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
