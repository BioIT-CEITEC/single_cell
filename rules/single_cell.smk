
# ###################################
# CELLRANGER RULES
#
import os

rule cellranger_call:
  input:
         fastq1=BR.remote(expand("raw_fastq/{lib}_{num}_R1.fastq.gz",zip,lib=LIBS,num=NUMS)),
         fastq2=BR.remote(expand("raw_fastq/{lib}_{num}_R2.fastq.gz",zip,lib=LIBS,num=NUMS)),
         cellranger_dir_tar = BR.remote(os.path.join(config["globalResources"], "tools/cellranger/cellranger-5.0.1.tar.gz")),
         transcriptome_files = BR.remote_input_dir(os.path.join(reference_directory, "tool_data/cellranger/refdata-gex-GRCh38-2020-A")),
         feature_ref_path = BR.remote(os.path.join(config["globalResources"], "tools/cellranger/feature_ref_files/feature_ref_" + f"{config['sc_hashtags']}.csv")),

  params:
          task_wdir = BR.get_path(os.path.join(config["globalTaskPath"], config["task_name"])),
          libraries = BR.get_path(os.path.join(config["entity_name"]+".lib.csv")),
          sample = SAMPLES,
          transcriptome_files_path = BR.get_path(os.path.join(reference_directory, "/tool_data/cellranger/refdata-gex-GRCh38-2020-A")),
          libs = list(library_types_dict.keys()),
          outdir = "cell_ranger",
          sc_hashtags = config["sc_hashtags"],
          library_types_dict = library_types_dict,
  output:
          directory(BR.get_path("cell_ranger/outs")),
          c1=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz",zip,lib=LIBS,num=NUMS)),
          c2=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz",zip,lib=LIBS,num=NUMS))

  log:    BR.remote("logs/all_samples/cellranger_call.log"),
  threads: 40
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"
