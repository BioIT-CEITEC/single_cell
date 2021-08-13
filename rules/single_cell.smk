
# ###################################
# CELLRANGER RULES
#

rule cellranger_call:
  input: c1 = expand("results/singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
         c2 = expand("results/singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS)
  params: sample = SAMPLES,
          libraries = os.path.join(config["library_name"]+".lib.csv"),
          wdir = os.getcwd(),
          outdir = "results/cell_ranger",
          binary = os.path.join(GLOBAL_REF_PATH,"general/cellranger/cellranger-5.0.1","cellranger"),
          feature_ref_dir = os.path.join(GLOBAL_REF_PATH,"general/cellranger/feature_ref_files"),
          sc_hashtags = config["sc_hashtags"],
          transcriptome = expand("{ref_dir}/other/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref=config["reference"])[0],
          HTML = "results/cell_ranger/outs/web_summary.html"
  output: REPORT = os.path.join("results",config["library_name"]+".final_report.html")
  log:   "logs/all_samples/cellranger_call.log",
  threads: 40
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"


####################################
# FASTQ_SYMLINK RULEs
#

rule fastq_symlink:
  input: fastq= "raw_fastq/{lib}_{num}{read_pair_tag}.fastq.gz",
  output: singleCell = "results/singleCell_fastq/{lib}/{lib}_S{num}_L001{read_pair_tag}_001.fastq.gz",
  log:    "logs/{lib}/{lib}_{num}{read_pair_tag}_singleCell_preprocess.log"
  threads: 1
  conda:  "../wrappers/fastq_symlink/env.yaml"
  script: "../wrappers/fastq_symlink/script.py"
