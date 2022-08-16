
# ###################################
# CELLRANGER RULES
#
import os

rule cellranger_call:
  input:
         fastq1=BR.remote(expand("raw_fastq/{lib}_{num}_R1.fastq.gz",zip,lib=LIBS,num=NUMS)),
         fastq2=BR.remote(expand("raw_fastq/{lib}_{num}_R2.fastq.gz",zip,lib=LIBS,num=NUMS)),
         cellranger_dir_zip = BR.remote(os.path.join(config["globalResources"], "tools/cellranger/cellranger-5.0.1.zip")),
         transcriptome_files = BR.remote(['/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/reference.json','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/pickle/genes.pickle','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/chrName.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/exonInfo.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/SAindex','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/Genome','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/geneInfo.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/sjdbList.out.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/sjdbInfo.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/transcriptInfo.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/sjdbList.fromGTF.out.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/chrLength.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/genomeParameters.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/chrStart.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/SA','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/chrNameLength.txt','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/star/exonGeTrInfo.tab','/resources/organisms/homo_sapiens/GRCh38/tool_data/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf']),
         libraries = BR.remote(os.path.join(config["entity_name"]+".lib.csv")),
         feature_ref_path = BR.remote(expand("{ref_path}/tools/cellranger/feature_ref_files/feature_ref_{sc_hashtag}.csv",
             sc_hashtag=config["sc_hashtags"], ref_path=config["globalResources"])),
  params:
          wdir = BR.remote(os.path.join(config["globalTaskPath"], config["task_name"])),
          sample = SAMPLES,
          outdir = "cell_ranger",
          transcriptome = BR.remote(expand("{ref_dir}/tool_data/cellranger/refdata-gex-{ref}",ref_dir=reference_directory,ref="GRCh38-2020-A")),
          sc_hashtags = config["sc_hashtags"],
          library_types_dict = library_types_dict,
  output:
          BR.remote("cell_ranger/outs/web_summary.html"),
          c1=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R1_001.fastq.gz",zip,lib=LIBS,num=NUMS)),
          c2=BR.remote(expand("singleCell_fastq/{lib}/{lib}_S{num}_L001_R2_001.fastq.gz",zip,lib=LIBS,num=NUMS))
  log:    BR.remote("logs/all_samples/cellranger_call.log"),
  threads: 40
  conda:   "../wrappers/cellranger_call/env.yaml"
  script:  "../wrappers/cellranger_call/script.py"
  # script: "../wrappers/cellranger_call/test_script.py"


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
