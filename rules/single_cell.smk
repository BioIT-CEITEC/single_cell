
# ###################################
# CELLRANGER RULES
#

#onerror:
#    print("Cellranger failed cause it's stupid!")
#    shell("ls -ld /tmp/*cell*")
#    shell("cat /tmp/cell_ranger/_log")

# rule cellranger_call:
#   input: c1 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
#          c2 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
#          tar = os.path.join(GLOBAL_REF_PATH,"general/cellranger/cellranger-"+config['cellranger_version']+".tar.gz")
#   output: "cell_ranger/outs/web_summary.html"
#   params: sample = SAMPLES,
#           libraries = config["entity_name"]+".lib.csv",
#           wdir = GLOBAL_TMPD_PATH,
#           outdir = "cell_ranger",
#           bin_path = os.path.join(GLOBAL_TMPD_PATH,"cellranger-"+config['cellranger_version']),
#           feature_ref_dir = os.path.join(GLOBAL_REF_PATH,"general/cellranger/feature_ref_files"),
#           sc_hashtags = config["sc_hashtags"],
#           library_types_dict = library_types_dict,
#           transcriptome = expand("{ref_dir}/other/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref=config["reference"])[0]
#   log:   "logs/all_samples/cellranger_call.log",
#   threads: 40
#   resources: mem = 100
#   conda:   "../wrappers/cellranger_call/env.yaml"
#   script:  "../wrappers/cellranger_call/script.py"

# ###################################
# STARsolo RULES
#

rule STARSolo_call:
    input:  c1 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
            c2 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
            genome = config["organism_fasta"],# defined in utilities
            gtf=config["organism_gtf"],# defined in utilities
            index=config["organism_star"]  # defined in utilities
    output: bam = "mapped/{sample}.solo.bam",
            bai = "mapped/{sample}.solo.bam.bai",
            # transcriptome_bam = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam",
            # transcriptome_bai = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam.bai",
    log:    "logs/{sample}/starsolo.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "mapped/{sample}/{sample}",
            strandness = config["strandness"],
            num_mismatch= 999,  # Maximum number of mismatches; set this to high number (999) to disable and to use only perc_mismatch
            perc_mismatch= config["perc_mismatch"],
            max_intron= config["max_intron"],# Default used by ENCODE is 1000000; to turn this off set 1
            max_mate_dist=1000000,# Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
            read_len=100,# Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
            organism=config["organism"],
            map_perc= config["map_perc"],
            map_score=config["map_score"],
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH,
    conda: "../wrappers/alignment_RNA/env.yaml"
    script: "../wrappers/alignment_RNA/script.py"


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
