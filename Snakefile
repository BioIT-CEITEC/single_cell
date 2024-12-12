import os
import json
import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoots utilities #####

module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

tools = BR.load_tooldir()

SAMPLES = [x for x in sample_tab.sample_name]
LIBS = [x.rsplit("_",1)[0] for x in SAMPLES]
NUMS = [x.rsplit("_",1)[1] for x in SAMPLES]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    lib_name="[^\.\/]+",
    read_pair_tag = "(_R.)?"

rule all:
    input:
        expand("mapped/{lib}_S{num}/{lib}_S{num}.rds", zip, lib = LIBS, num = NUMS, sample = sample_tab.sample_name),
        expand("mapped/{lib}_S{num}/{lib}_S{num}_qc_report.html", zip, lib = LIBS, num = NUMS, sample = sample_tab.sample_name)

##### Modules #####

include: "rules/single_cell.smk"
