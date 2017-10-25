# velocyto.R
RNA velocity estimation in R

## Installation
The easiest way to install velocyto.R is using devtools::install_github() from R:
```
library(devtools)
install_github("velocyto-team/velocyto.R")
```

## Tutorials

### [Chromaffin / SMART-seq2](http://pklab.med.harvard.edu/velocyto/notebooks/R/chromaffin.nb.html)
The example shows how to annotate SMART-seq2 reads from bam file and estimate RNA velocity.

### [Dentate Gyrus / loom](http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html)
The example shows how to load spliced/unspliced matrices from loom files prepared by [velocyto.py CLI](http://velocyto.org/velocyto.py/tutorial/index.html#running-the-cli), use [pagoda2](https://github.com/hms-dbmi/pagoda2) to cluster/embed cells, and then visualize RNA velocity on that embedding.

### [Mouse BM / dropEst](http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html)
This example shows how to start analysis using dropEst count matrices, which can calculated from inDrop or 10x bam files using [dropEst pipeline](https://github.com/hms-dbmi/dropEst/). It then uses [pagoda2](https://github.com/hms-dbmi/pagoda2) to cluster/embed cells, and then visualize RNA velocity on that embedding.
