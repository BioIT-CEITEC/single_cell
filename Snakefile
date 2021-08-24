import os
import json
import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"


if not "sc_hashtags" in config:
    config["sc_hashtags"] = "no"


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]


##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

SAMPLES = [x for x in sample_tab.sample_name]
LIBS = [x.rsplit("_",1)[0] for x in SAMPLES]
NUMS = [x.rsplit("_",1)[1] for x in SAMPLES]

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
else:
    read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    lib_name="[^\.\/]+",
    read_pair_tag = "(_R.)?"



rule all:
    input: os.path.join(config["library_name"]+".final_report.html")

##### Modules #####

include: "rules/single_cell.smk"




