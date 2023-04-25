# Papet

---------------------------------------------

This repository contains the PAcbio Prediction of Epigenetics Technology (Papet) software to model the PacBio CCS kinetic signal and predict epigenetics modifications from this kinetic signal.


## Table of content

1. [Content](#content)
2. [Dependencies](#dependencies)
3. [Docker image](#docker-image)
4. [Compilation](#compilation)
5. [About PacBio kinetics and epigenetics](#about-pacbio-kinetics-and-epigenetics)
6. [Running Papet](#running-papet)  
    6.1. [model-kinetic](#kinetics)  
    6.2. [model-kinetic-txt](#model-kinetic-txt)  
    6.3. [kinetics](#kinetics)  
    6.4. [kinetics-wig](#kinetics-wig)  
    6.5. [kinetics-kmer](#kinetics-kmer)  
    6.6. [model-sequence](#model-sequence)  
    6.7. [model-sequence-txt](#model-sequence-txt)  
    6.8. [predict](#predict)
7. [Acknowledgments](#acknowledgments)

## Dependencies

Papet relies on the following third party libraries:

- pthread, normally a system library on linux system.
- [boost](https://www.boost.org/) v1.78 or higher
- [pbcopper](https://github.com/PacificBiosciences/pbcopper) v2.0.0 or higher
- [pbbam](https://github.com/PacificBiosciences/pbbam) v2.0.0 (pbbam relies on pbcopper)
- [ngsaipp](https://github.com/ngs-ai-org/ngsaipp)

We strongly recommend to install these libraries in `/usr/local/lib` and their corresponding header files in `/usr/local/include`. If done this way, cmake compilation configuration should work out of the box.

## Docker image

A Dockerfile is present in the repository to build a docker image. To do so, simply use:
```
docker build -t papet .
```

Alternatively, the already built image can be pulled from [dockerhub](https://hub.docker.com/r/ngsai/papet) using
```
docker pull ngsai/papet
```

## Compilation

The building process uses cmake. To compile and install papet, simply type:
```
./build.sh
```

The `CMakeLists.txt` file contains 3 variable defined at its top:

```
# user defined paths
## list of directories containing libraries headers
set(INCLUDE_DIRECTORIES "/usr/local/include/"
                        "/usr/local/include/ngsaipp")
## list of directories in which the required libraries are installed
set(LINK_DIRECTORY      "/usr/local/lib")
## path in which the final executable will be installed
set(INSTALL_DIRECTORY   "/usr/local/bin")
```

- `INCLUDE_DIRECTORIES` a space separated list of directories in which the required third party library header files are located.

- `LINK_DIRECTORY` a space separated list of directories in which all the required libraries are installed.

- `INSTALL_DIRECTORY` the directory in which papet will be installed. 


## About PacBio kinetics and epigenetics

PacBio BAM file epignetics signal data handling,  modelling as well as the epigenetics prediction computations are performed by ngsaipp. For an in-depth dive on this topic, please read [the dedicated documention](https://github.com/ngs-ai-org/ngsaipp/blob/master/PacBio_kinetics.md)


## Running Papet

Papet is a wrapper application that can run different commands to fullfil specific tasks. The general synthax is:

```
papet [command] [options]
```

Papet accepts the following options and commands:

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-vaudois           | Tastes quite good.|
  |       | model-kinetic         | Creates kinetic signal models from CCSs. |
  |       | model-kinetic-txt     | Dumps a kinetic signal model in txt format. |
  |       | kinetics              | Extracts CCS kinetic information in txt format. |
  |       | kinetics-wig          | Creates WIG tracks from CCSs. |
  |       | kinetics-kmer         | Computes the per-kmer distribution of kinetic signal from CCSs. |
  |       | model-sequence        | Creates DNA sequence kinetic signal models from CCSs. |
  |       | model-sequence-txt    | Dumps a DNA sequence kinetic signal model in txt format. |
  |       | predict               | Predicts the presence of epignetic modifications from CCSs. |


### model-kinetic

model-kinetic trains a kinetic signal model by computing the kinetic (IPD and PWD) signal distribution at each position inside windows of \<N\> bases. The signal distributions are learned by extracting the kinetic signal, from the reads, using \<N\> bases long windows centered on the given CpGs. The model is then serialized in the given file. The final results is a serialized instance derived from the KineticModel class ([see ngsaipp KineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/KineticModel.hpp)).

The synthax is:

```
papet model-kinetic [type] [options]
```

The exact type of kinetic signal model is defined using the first argument. The accepted values are:

- `raw`: simply computes the distribution of IPD and PWD at each position in the window. The model is a serialized instance of the RawKineticModel class ([see ngsaipp RawKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/RawKineticModel.hpp)).
- `raw-norm`: computes the normalized IPD and PWD' distributions at each position in the window. The model is a serialized instance of the NormalizedKineticModel class ([see ngsaipp NormalizedKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/NormalizedKineticModel.hpp)).
- `diposition`: computes the distributions of IPD and PWD signal at each position in the window as a function of its direct neighbour. The model is a serialized instance of the DiPositionKineticModel class ([see ngsaipp DiPositionKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/DiPositionKineticModel.hpp)).
- `diposition-norm`: computes the distributions of normalized IPD and PWD signal at each position in the window as a function of its direct neighbor. The model is a serialized instance of the DiPositionNormalizedKineticModel 
class ([see ngsaipp DiPositionNormalizedKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp)).
- `pairwise`: computes all distributions of IPD and PWD as a function of the signal at another position in the window. The model is a serialized instance of the PairWiseKineticModel class ([see ngsaipp PairWiseKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/PairWiseKineticModel.hpp))
- `pairwise-norm`: computes all distributions of normalized IPD and PWD as a function of the signal at another position in the window. The model is a serialized instance of the PairWiseNormalizedKineticModel class ([see ngsaipp PairWiseNormalizedKineticModel.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp))


This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam               | A coma separated list of paths to the bam files containing the mapped PacBio CCS of interest. |
  |       | \-\-bed               | The path to the bed file containing the genomic coordinates of the CpGs of interest.|
  |       | \-\-out               | The path to file in which the kinetic model will be saved.|
  |       | \-\-background       | For normalized models only, the path to background model to use. It must contain a serialized KmerMap.|
  |       | \-\-size              | The size of the model, the length of the signal window to model, in bp|
  |       | \-\-nbin              | The number of bins in each histogram. |
  |       | \-\-xmin              | The lower limit of the lower bin in each histogram. |
  |       | \-\-xmax              | The upper limit of the upper bin in each histogram.|
  |       | \-\-pseudocount       | A number of counts that will be added to each bin in each histogram, by default 0.|
  |       | \-\-thread            | The number of threads, by default 1. |


### model-kinetic-txt

model-kinetic-txt converts a kinetic signal model file in tsv format and prints the results on stdout. If the model is a normalized one, the KmerMap is not tincluded in the tsv conversion.

The synthax is:

```
papet model-kinetic-txt [options] [> FILE]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-model             | The path to the file containing the kinetic  model to convert in tsv format. |


### kinetics

kinetics is an application to extract interpulse duration (IPDs) and pulse widths (PWs) from mapped PacBio CCS reads that overlap a given set of genomic regions specified in a BED 6 file. The results are printed on stdout in tsv format. The first row is a header. Then, each line contains per read 
sequence, IPDs and PWDs.
The memory footprint of this progam is **O(n)**, where **n** is the window size. The results are immediately printed and nothing is kept in memory.

The synthax is:

```
papet kinetics [options] [> FILE]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam   arg         | A coma separated list of paths to the bam files containing the mapped PacBio CCS of interest.|
  |       | \-\-bed  arg          | The path to the bed file containing the genomic regions of interest. |
  |       | \-\-winSize           | The size in bp of the windows from which the CCS features will be extracted. It must be odd. The window will be centered on the center of the genomic regions of interst. |


### kinetics-wig

kinetics-wig is an application that creates wig kinetic tracks from a set of PacBio CCS from from bam files, over a set of defined genomic windows. Exactly 4 tracks are created each time, corresponding to the IPD on forward strand, the IPD on the reverse strand, the PWD on the forward strand and the PWD on the reverse strand.

The synthax is:

```
papet kinetics-wig [options]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam arg         | A coma separated list of paths to the bam files containing the mapped PacBio CCS from which the IPD/PWD track must be computed.|
  |       | \-\-bed               |  The path to the bed file containing the genomic oordinates of the CpGs of interest. |
  |       | \-\-out               |  A path prefix to use to write the results files. In total, 4 resulting files will be created with this prefix : <prefix>_IPDfw.wig, <prefix>_IPDrv.wig, <prefix>_PWDfw.wig and <prefix>_PWDrv.wig containing the IPD and PWD forward and reverse track respectively. |
  |       | \-\-winSize           |  The size of the window (in bp) around the CpGs in which the average kinetic signal will be computed. |


### kinetics-kmer

kinetics-kmer is an application that computes the interpulse duration (IPDs) and pulse widths (PWs) value distribution for all kmers found in a PacBio CCS dataset. The kmer IPD / PWD value is defined has the IPD / PWD value  corresponding to the central position of the kmer.
The results are printed on stdout in tsv format. The first row is a header. Then, each line contains the results for one kmer. Each line contains the kmer sequence (column 1), the IPD distribution for this kmer (columns 2 to 954) amd the PWD distribution (columns 955 1907). IPD / PWD values can take any score in the [0,952] interval. Thus there are 953 columns for each IPD / PWD distribution. The 1st IPD / PWD related value (column 2 / 955) corresponds to the number of times an IPD / PWD value of 0 was observed for this kmer. The 2nd IPD / PWD value (column 3 / 956) corresponds to the number of time an IPD / PWD value of 1 was observed for this kmer, and so on until a score of 952 (column 953 / 1907). In the final table, not all possible kmers will be present in the table, only the ones found in
the data. The memory footprint of this progam is **O(4^k)** where **k** is the size of the kmers.

The synthax is:

```
papet kinetics-kmer [options]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam   arg         | A coma separated list of paths to the bam files containing the PacBio CCS of interest.|
  |       | \-\-kmer arg          | The size of the kmer. It must be odd. |


### model-sequence

model-sequence is an application that computes the mean IPD and PWD per position signal for all possible kmer from a set of CCS reads (see [Expected kinetic background](https://github.com/ngs-ai-org/ngsaipp/blob/master/PacBio_kinetics.md#background-kinetic-model)). The final results is a 
serialized instance of KmerMap ([see ngsaipp KmerMap.hpp](https://github.com/ngs-ai-org/ngsaipp/tree/master/include/ngsaipp/epigenetics/KmerMap.hpp)). The serialized instance can then either be loaded in another program or converted in tsv format using `model-sequence-txt`.
The memory footprint of this progam is **O(4^k)** where 
**k** is the size of the kmers. Note that overflow errors may occure if a given kmer is encountered too many times in the dataset. This can be the 
case if there are too many reads (unlikely) or if the kmer size is too small (experience shows that kmer of 3 leads to this issue). In any case, if the mean signal seems weird, this is a likely explanation.

The synthax is:
```
papet model-sequence [options]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam   arg         | A coma separated list of paths to the bam files containing the PacBio CCS of interest.|
  |       | \-\-out   arg         | The path to the file in which the KmerMap will be dumped.|
  |       | \-\-kmer arg          | The size of the kmer. It must be odd. |


### model-sequence-txt

model-sequence-txt is an application that converts into tsv format a serialized KmerMap created by `model-sequence`. The results are printed on stdout. The first row contains headers, the following rows contains the kmer data. Each row contains the kmer sequence, the kmer number of occurences, `k` per position IPD means and `k` per position PWD means. 
The memory footprint of this progam is **O(4^k)** where **k** is the size of the kmers. Note that overflow errors may occure if a given kmer is encountered too many times in the dataset. This can be the case if there are too many reads (unlikely) or if the kmer size is too small (experience shows that kmer of 3 leads to this issue). In any case, if the mean signal seems weird, this is a likely explanation.

The synthax is:
```
papet model-sequence-txt [options] [>FILE]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-model arg         | The path to the file containing the KmerMap to convert in tsv format.|


### predict

predict predicts the presence of epigenetic modifications at the CpG of interest given two models and returns the results on stdout in BED 6 format. The score field contains the probability of the presence of an epigenetic modification.

The synthax is:
```
papet predict [options] [>FILE]
```

This program has the following options :

  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:----------------------|:--------------------------|
  | -h    | \-\-help              | Produces the help message |
  |       | \-\-bam               | A coma separated list of paths to the bam files containing the mapped PacBio CCS of interest. |
  |       | \-\-bed               | The path to a bed file containing the coordinates of the CpGs interest. |
  |       | \-\-modelMeth         | The path to the file containing the methylated kinetic model to use.  |
  |       | \-\-modelUnmeth       | The path to the file containing the unmethylated kinetic model to use.  |
  |       | \-\-prob              | The prior probability of methylation for any CpG. It must belong to [0,1]. 0.5 by default.  |
  |       | \-\-thread            | The number of threads, by default 1.  |


## Acknowledgments

Richard Hall and Pacific Biosciences for the discussions and the shared data.

## Authors

* **Romain Groux**
