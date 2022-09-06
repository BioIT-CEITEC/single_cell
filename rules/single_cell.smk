
# ###################################
# CELLRANGER RULES
#
import os

rule cellranger_call:
  input:
         fastq1=BR.remote(expand("raw_fastq/{lib}_{num}_R1.fastq.gz",zip,lib=LIBS,num=NUMS)),
         fastq2=BR.remote(expand("raw_fastq/{lib}_{num}_R2.fastq.gz",zip,lib=LIBS,num=NUMS)),
         cellranger_dir_zip = BR.remote(os.path.join(config["globalResources"], "tools/cellranger/cellranger-5.0.1.zip")),
         transcriptome_files = BR.remote_input_dir(expand("{ref_dir}/tool_data/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref="GRCh38-2020-A")),
         feature_ref_path = BR.remote(expand("{ref_path}/tools/cellranger/feature_ref_files/feature_ref_{sc_hashtag}.csv",
             sc_hashtag=config["sc_hashtags"], ref_path=config["globalResources"])),
  params:
          wdir = BR.remote(os.path.join(config["globalTaskPath"], config["task_name"])),
          libraries = BR.remote(os.path.join(config["entity_name"]+".lib.csv")),
          sample = SAMPLES,
          transcriptome_path = BR.remote(expand("{ref_dir}/tool_data/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref="GRCh38-2020-A")),
          libs = list(library_types_dict.keys()),
          outdir = "cell_ranger",
          sc_hashtags = config["sc_hashtags"],
          library_types_dict = library_types_dict,
  output:
          directory(BR.remote("cell_ranger/outs")),
          c1=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz",zip,lib=LIBS,num=NUMS)),
          c2=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz",zip,lib=LIBS,num=NUMS))

  log:    BR.remote("logs/all_samples/cellranger_call.log"),
  threads: 40
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"


# rule test:
#       input: BR.remote_dir(os.path.join(config["globalTaskPath"],config["task_name"],"test"))
#       output: res = directory(BR.remote_dir("result_dir"))
#       script: "../wrappers/cellranger_call/test_script.py"


  ####################################
# FASTQ_SYMLINK RULEs
#
# rule fastq_symlink:
#   input: fastq= BR.remote("raw_fastq/{lib}_{num}{read_pair_tag}.fastq.gz"),
#   output: singleCell = BR.remote("singleCell_fastq/{lib}/{lib}_S{num}_L001{read_pair_tag}_001.fastq.gz"),
#   log:    BR.remote("logs/{lib}/{lib}_{num}{read_pair_tag}_singleCell_preprocess.log"),
#   params: wdir = BR.remote(os.path.join(config["globalTaskPath"], config["task_name"])),
#   threads: 1
#   conda:  "../wrappers/fastq_symlink/env.yaml"
#   script: "../wrappers/fastq_symlink/script.py"
