import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

if not "sc_hashtags" in config:
    config["sc_hashtags"] = "no"

if not 'cellranger_version' in config:
    config['cellranger_version'] = "5.0.1"

##### Config processing #####

# Samples
#

SAMPLES = [x for x in sample_tab.sample_name]
LIBS = [x.rsplit("_",1)[0] for x in SAMPLES]
library_types_dict = {}
for idx, val in enumerate(sample_tab.SC_lib_type):
    if LIBS[idx] in library_types_dict and library_types_dict[LIBS[idx]] != val:
        raise ValueError("All 'SC type' for one SC library must be the same (Gene Expression or Antibody Capture). Fix library setting.")
    else:
        library_types_dict[LIBS[idx]] = val
NUMS = [x.rsplit("_",1)[1] for x in SAMPLES]


wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    lib_name="[^\.\/]+",
    read_pair_tag = "(_R.)?"



rule all:
    input: "cell_ranger/outs/web_summary.html"

##### Modules #####

include: "rules/single_cell.smk"




