
# ###################################
# CELLRANGER RULES
#

#onerror:
#    print("Cellranger failed cause it's stupid!")
#    shell("ls -ld /tmp/*cell*")
#    shell("cat /tmp/cell_ranger/_log")

rule cellranger_call:
  input: c1 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
         c2 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
         tar = os.path.join(GLOBAL_REF_PATH,"general/cellranger/cellranger-"+config['cellranger_version']+".tar.gz")
  output: "cell_ranger/outs/web_summary.html"
  params: sample = SAMPLES,
          libraries = config["entity_name"]+".lib.csv",
          wdir = GLOBAL_TMPD_PATH,
          outdir = "cell_ranger",
          bin_path = os.path.join(GLOBAL_TMPD_PATH,"cellranger-"+config['cellranger_version']),
          feature_ref_dir = os.path.join(GLOBAL_REF_PATH,"general/cellranger/feature_ref_files"),
          sc_hashtags = config["sc_hashtags"],
          library_types_dict = library_types_dict,
          transcriptome = expand("{ref_dir}/other/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref=config["reference"])[0]
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
