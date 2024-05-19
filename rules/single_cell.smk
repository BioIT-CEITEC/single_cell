
# ###################################
# CELLRANGER RULES
#
rule cellranger_call:
  input: c1 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
         c2 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
         tar = config["tool_dir"] + "/cellranger/cellranger-" + config['cellranger_version']+ ".tar.gz"
  output: "cell_ranger/outs/web_summary.html"
  params: sample = SAMPLES,
          libraries = config["entity_name"]+".lib.csv",
          wdir = GLOBAL_TMPD_PATH,
          outdir = "cell_ranger",
          bin_path = os.path.join(GLOBAL_TMPD_PATH,"cellranger-"+config['cellranger_version']),
          feature_ref_dir = config["tool_dir"] + "/cellranger/feature_ref_files",
          sc_hashtags = config["sc_hashtags"],
          library_types_dict = library_types_dict,
          transcriptome = config["organism_transcriptome"], #defined in bioroots utilities
  log:   "logs/all_samples/cellranger_call.log",
  threads: 40
  resources: mem = 100
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"


####################################
# FASTQ_SYMLINK RULEs
#

rule fastq_symlink:
  input: fastq= "raw_fastq/{lib}_{num}{read_pair_tag}.fastq.gz",
  output: singleCell = "singleCell_fastq/{lib}/{lib}_S{num}_L001{read_pair_tag}_001.fastq.gz",
  log:    "logs/{lib}/{lib}_{num}{read_pair_tag}_singleCell_preprocess.log",
  params: wdir = os.getcwd()
  threads: 1
  conda:  "../wrappers/fastq_symlink/env.yaml"
  script: "../wrappers/fastq_symlink/script.py"
