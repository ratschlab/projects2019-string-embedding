# projects2019-string-embedding
This project contains codes that simulate simple whole-genome evolution, and several analysis methods to find the MSA blocks.

# requirements to run the python codes
The python code must be run by Python 3.6+. The packages required are `shutil`, `annoy`, `tqdm`, and `flann`, which all can be installed by conda. 

For the c++ codes, they have to be compiled with `c++14`, but preferably `c++17` in case of future code written with the updated standards. 

# config file

put a folder named ".config_string_embedding" in your home  directory that contains the path to the project data directory under the variable `PROJ_DIR`. Currently the data is located on the cluster at "/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic" and so the config file should look like:

`PROJ_DIR = /cluster/work/grlab/projects/projecs2019-string-embedding/synthetic`

if you would like to run the codes locally, it `PROJ_DIR` must be equal to the directory that contains the sample `/data` directory in the project. These datasets are small and only for testing purposes. 

Currently, under `$(PROJ_DIR)` there are two main directories `$(PROJ_DIR)/data$` and `$(PROJ_DIR)/results`. Unless specified by a full path, the results will be saved in the results directory
