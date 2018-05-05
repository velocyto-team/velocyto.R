# velocyto.R
RNA velocity estimation in R

## System requirements
velocyto.R can be installed on unix-flavored systems, and requires the following key elements:

* C++11
* Open MP support
* boost libaries
* igraph library
* hdf5c++ library (as required by the h5 R package to support loom files)

## Installation
The easiest way to install velocyto.R is using devtools::install_github() from R:
```
library(devtools)
install_github("velocyto-team/velocyto.R")
```
You need to have boost (e.g. `sudo apt-get install libboost-dev`) and openmp libraries installed. You can see detailed installation commands in the dockers/debian9/Dockerfile. 

### Dockers
If you are having trouble installing the package on your system, you can build a docker instance that can be used on a wide range of systems and cloud environments. To install docker framework on your system see [installation instruction](https://github.com/wsargent/docker-cheat-sheet#installation). After installing the docker system, use the following commands to build a velocyto.R docker instance:
```bash
cd velocyto.R/dockers/debian9
docker build -t velocyto .
docker run --name velocyto -it velocyto
```

## Tutorials

### [Chromaffin / SMART-seq2](http://pklab.med.harvard.edu/velocyto/notebooks/R/chromaffin2.nb.html)
The example shows how to annotate SMART-seq2 reads from bam file and estimate RNA velocity.

### [Dentate Gyrus / loom](http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html)
The example shows how to load spliced/unspliced matrices from loom files prepared by [velocyto.py CLI](http://velocyto.org/velocyto.py/tutorial/index.html#running-the-cli), use [pagoda2](https://github.com/hms-dbmi/pagoda2) to cluster/embed cells, and then visualize RNA velocity on that embedding.

### [Mouse BM / dropEst](http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html)
This example shows how to start analysis using dropEst count matrices, which can calculated from inDrop or 10x bam files using [dropEst pipeline](https://github.com/hms-dbmi/dropEst/). It then uses [pagoda2](https://github.com/hms-dbmi/pagoda2) to cluster/embed cells, and then visualize RNA velocity on that embedding.
