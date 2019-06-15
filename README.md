# projects2019-string-embedding
This project contains codes that simulate simple whole-genome evolution, and several analysis methods to find the MSA blocks. 

# config file

put a folder named ".config_string_embedding" in your home  directory that contains the path to the project data directory under the variable `PROJ_DIR`. Currently the data is located at "/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic" and so the config file should look like:

`PROJ_DIR = /cluster/work/grlab/projects/projecs2019-string-embedding/synthetic`

Currently, under `$(PROJ_DIR)` there are two main directories `$(PROJ_DIR)/data$` and `$(PROJ_DIR)/results`. Unless specified by a full path, the results will be saved in the results directory
