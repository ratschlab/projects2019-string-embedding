# projects2019-string-embedding
This project implements a fast sequence similarity sketching, with experiments to compare it against min-hash and ordered min-hash. Here is a priliminary result comparing this method (upper plots) to ordered min-hash (lower plots):
![The main result](https://github.com/ratschlab/projects2019-string-embedding/blob/master/result.png)

## organization
The overal organization of the code is as follows:
- `$(HOME)/.config_string_embedding`: config file of the project, specifying `PROJECT_DIR`
- code: The python codes are mostly for prototyping and eventually all the codes must run natively in C++ (in progress).
  - python: *.py files for sequence genration, string embedding, and and simple KNN search
  - cpp: *.cpp files for sequence sketching  (in progress)
  - matlab: *.m containing the lastest development 
- data: sample data for unit testing
- `$(PROJ_DIR)`
  - data
  - results 


## requirements to run the codes
The python code must be run by Python 3.6+. The packages required are `numpy, shutil, annoy, tqdm`, and `flann`, which all can be installed by conda.  

For the c++ codes, the code has been compiled and tested with `-std=c++17` standards. There is cmake to compile the cpp files automatically.  


## config file

put a folder named `.config_string_embedding` in your home  directory that contains the path to the project data directory under the variable `PROJ_DIR`. Currently the data is located on the cluster at `/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic` and so the config file should look like

`PROJ_DIR = /cluster/work/grlab/projects/projecs2019-string-embedding/synthetic`

if you would like to run the codes locally, it `PROJ_DIR` must be equal to the directory that contains the sample `/data` directory in the project. These datasets are small and only for testing purposes. 

Currently, under `$(PROJ_DIR)` there are two main directories `$(PROJ_DIR)/data$` and `$(PROJ_DIR)/results`. Unless specified by a full path, the results will be saved in the results directory
