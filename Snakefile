from snakemake.utils import min_version

min_version("7.2.1")

configfile: "config.json"

config["computing_type"] = "kubernetes"

module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="kube_dirs")
    config: config

use rule * from BR as other_*

#### Config preprocessing ####

if not "sc_hashtags" in config:
    config["sc_hashtags"] = "no"

sample_tab = BR.load_sample()
read_pair_tags = BR.set_read_pair_tags()

config["lib_ROI"] = "RNA"

#### Reference processing ####

config = BR.load_ref()
config = BR.load_organism()
reference_directory = BR.reference_directory()

#### Samples processing ####

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
    input: BR.remote(expand("cell_ranger/outs")),
    output: BR.remote("completed.txt")
    shell: "touch {output}"

##### Modules #####

include: "rules/single_cell.smk"




