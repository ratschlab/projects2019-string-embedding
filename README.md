# projects2019-string-embedding
This project contains codes that simulate simple whole-genome evolution, and several analysis methods to find the MSA blocks.

The overal organization of the code is as follows:
- `$(HOME)/.config_string_embedding`: config file of the project, specifying `PROJECT_DIR`
- code
  - python: *.py files for sequence genration, string embedding, and KNN search
  - cpp: *.cpp files for sequence genration and Ordered-LSH (in progress)
- data: sample data for unit testing
- `$(PROJ_DIR)`
  - data
  - results 


# requirements to run the python codes
The python code must be run by Python 3.6+. The packages required are `numpy, shutil, annoy, tqdm`, and `flann`, which all can be installed by conda.  

For the c++ codes, the code has been compiled and tested with `-std=c++17` standards. 

The python codes are mostly for prototyping and eventually all the codes must run natively in C++ (in progress).

# config file

put a folder named `.config_string_embedding` in your home  directory that contains the path to the project data directory under the variable `PROJ_DIR`. Currently the data is located on the cluster at `/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic` and so the config file should look like

`PROJ_DIR = /cluster/work/grlab/projects/projecs2019-string-embedding/synthetic`

if you would like to run the codes locally, it `PROJ_DIR` must be equal to the directory that contains the sample `/data` directory in the project. These datasets are small and only for testing purposes. 

Currently, under `$(PROJ_DIR)` there are two main directories `$(PROJ_DIR)/data$` and `$(PROJ_DIR)/results`. Unless specified by a full path, the results will be saved in the results directory
