import os
import json
import pandas as pd
from snakemake.utils import min_version

min_version("7.2.1")

configfile: "config.json"

config["computing_type"] = "kubernetes"
GLOBAL_REF_PATH = "/mnt/references/"

module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="main")
    config: config

use rule * from BR as other_*

if not "sc_hashtags" in config:
    config["sc_hashtags"] = "no"

# setting organism from reference
#f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
#reference_dict = json.load(f)
#f.close()
#config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
sample_tab = BR.load_sample()
config["lib_ROI"] = "RNA"
BR.load_ref()
BR.load_organism()


##### Config processing #####
#
# Folders
#
reference_directory = BR.reference_directory()

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

# if not config["is_paired"]:
#     read_pair_tags = [""]
# else:
#     read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    lib_name="[^\.\/]+",
    read_pair_tag = "(_R.)?"



rule all:
    input: BR.remote("cell_ranger/outs/web_summary.html")

##### Modules #####

include: "rules/single_cell.smk"




