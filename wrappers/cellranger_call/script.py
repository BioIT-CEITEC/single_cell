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

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cellranger_call \n##\n")
f.close()


###### FASTQ preprocessing
# cellranger requires libraries to conform to a specific naming scheme

# create target directories
for singlecell_fastq in set(map(os.path.dirname, snakemake.output.c1 + snakemake.output.c2)):
    command = "mkdir -p " + singlecell_fastq
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)


# copy and rename fastq
def copy_fastq(raw_fastqs, output_fastqs):
    assert len(raw_fastqs) == len(output_fastqs)
    for i in range(len(raw_fastqs)):
        command = "cp " + raw_fastqs[i] + " " + output_fastqs[i]
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)


copy_fastq(snakemake.input.fastq1, snakemake.output.c1)
copy_fastq(snakemake.input.fastq2, snakemake.output.c2)


# CREATE the csv file if not existing and add the header
lf = open(snakemake.input.libraries, "w")
lf.write("fastqs,sample,library_type\n")
lf.close()
f = open(log_filename, 'at')
f.write("## COMMAND: create" + snakemake.input.libraries + " file\n")
f.close()

# ADD info to csv
# CHECK if "GE" is in the name of the library
for x in snakemake.params.libs:
    f = open(log_filename, 'at')
    f.write("## COMMAND: create " + x + " folder\n")
    f.close()

    SAMPLE_LIB_DIR = os.path.join("/tmp", snakemake.params.wdir, "singleCell_fastq", x)
    line_to_write = SAMPLE_LIB_DIR + "," + x + "," + snakemake.params.library_types_dict[x] + "\n"
    lf = open(snakemake.input.libraries, "at")
    lf.write(line_to_write)
    lf.close()

f = open(log_filename, 'at')
f.write("## COMMAND: filling" + snakemake.input.libraries + " file\n")
f.close()

if snakemake.params.sc_hashtags != "no":
    feature_ref_parameter = "--feature-ref=" + "/tmp/" + snakemake.input.feature_ref_path
else:
    feature_ref_parameter = ""

# CALL cellranger
cmd = "chmod +x " + snakemake.params.wdir + "/cellranger-5.0.1/bin/cellranger"
shell(cmd)

command = "cd " + snakemake.params.wdir + "; rm -Rf " + snakemake.params.outdir + " ; " + "cellranger-5.0.1/bin/cellranger" + " count " + \
          " --id=" + snakemake.params.outdir + \
          " --libraries=" + os.path.basename(snakemake.input.libraries) + \
          " " + feature_ref_parameter + \
          " --transcriptome=" + "/tmp/" + snakemake.params.transcriptome_path + \
          " --localcores " + str(snakemake.threads) + \
          " >> " + log_filename.replace(snakemake.params.wdir + "/", "") + " 2>&1 ; cd .."


f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
