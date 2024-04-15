rule STARSolo_call:
    input:  c1 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
            c2 = expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz", zip,  lib=LIBS, num=NUMS),
            genome = config["organism_fasta"],# defined in utilities
            gtf=config["organism_gtf"],# defined in utilities
            index=config["organism_starsolo"],  # defined in utilities
            whitelist = CR_whitelist
    output: bam = "mapped/{lib}_S{num}.solo.bam",
            bai = "mapped/{lib}_S{num}.solo.bam.bai",
    log:    "logs/{sample}/starsolo.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "mapped/{lib}_S{num}/{lib}_S{num}",
            tmpd = GLOBAL_TMPD_PATH,
            assay = config["assay_type"],
            tool_folder = config["tool_path"]
    conda: "../wrappers/alignment_RNA/env.yaml"
    script: "../wrappers/alignment_RNA/script.py"

rule fastq_symlink:
  input: fastq= "raw_fastq/{lib}_{num}{read_pair_tag}.fastq.gz",
  output: singleCell = "singleCell_fastq/{lib}/{lib}_S{num}_L001{read_pair_tag}_001.fastq.gz",
  log:    "logs/{lib}/{lib}_{num}{read_pair_tag}_singleCell_preprocess.log",
  params: wdir = os.getcwd()
  threads: 1
  conda:  "../wrappers/fastq_symlink/env.yaml"
  script: "../wrappers/fastq_symlink/script.py"

rule seurat_obj:
  input:  counts = "mapped/{lib}_S{num}/{lib}_S{num}_Solo.out/Gene/filtered/barcodes.tsv"
  output: rds = "mapped/{lib}_S{num}/{lib}_S{num}.rds",
          folder = "mapped/{lib}_S{num}/Plots/{lib}_S{num}_1_UMAP.png"
  log:   "logs/{lib}_S{num}/{lib}_S{num}_seurat.log",
  conda:   "../wrappers/seurat/env.yaml"
  script:  "../wrappers/seurat/script.py"

rule qc_report:
  input: csv = "mapped/{lib}_S{num}/{lib}_S{num}_Solo.out/Gene/Summary.csv"
  output: html = "mapped/{lib}_S{num}/{lib}_S{num}_qc_report.html"
  log:   "logs/{lib}_S{num}/{lib}_S{num}_qc_report.log"
  conda:  "../wrappers/qual_report/env.yaml"
  script:  "../wrappers/qual_report/script.py"