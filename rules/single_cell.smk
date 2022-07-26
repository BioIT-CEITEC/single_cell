
# ###################################
# CELLRANGER RULES
#
import os

rule cellranger_call:
  input: c1 = BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS)),
         c2 = BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS)),
         binary = BR.remote(os.path.join(GLOBAL_REF_PATH,"general/cellranger/cellranger-5.0.1","cellranger")),
         libraries = BR.remote(os.path.join(config["entity_name"]+".lib.csv")),
  params: sample = SAMPLES,
          outdir = "cell_ranger",
          wdir = BR.remote(os.path.join(config["globalTaskPath"], config["task_name"])),
          sc_hashtags = config["sc_hashtags"],
          transcriptome = BR.remote(expand("{ref_dir}/other/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref=config["reference"])[0]),
          feature_ref_dir = BR.remote(os.path.join(GLOBAL_REF_PATH,"general/cellranger/feature_ref_files")),
          library_types_dict = library_types_dict,
  output: BR.remote("cell_ranger/outs/web_summary.html"),
  log:   BR.remote("logs/all_samples/cellranger_call.log"),
  threads: 40
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"


####################################
# FASTQ_SYMLINK RULEs
#

rule fastq_symlink:
  input: fastq= BR.remote("raw_fastq/{lib}_{num}{read_pair_tag}.fastq.gz"),
  output: singleCell = BR.remote("singleCell_fastq/{lib}/{lib}_S{num}_L001{read_pair_tag}_001.fastq.gz"),
  log:    BR.remote("logs/{lib}/{lib}_{num}{read_pair_tag}_singleCell_preprocess.log"),
  params: wdir = os.getcwd()
  threads: 1
  conda:  "../wrappers/fastq_symlink/env.yaml"
  script: "../wrappers/fastq_symlink/script.py"
