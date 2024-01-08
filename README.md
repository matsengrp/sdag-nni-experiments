# NNI-Search Experiments

To view the results of the experiments, see the Jupyter notebooks 
`plot_generalized_pruning_results.ipynb`  and `plot_top_pruning_results.ipynb`.

## Dependencies and prerequisities

There are several dependencies and requirements before running an nni-search.

First, install bito: [bito](https://github.com/phylovi/bito/tree/main). 

Second, update your bito environment based on the [watching-mb repo](https://github.com/matsengrp/watching-mb).
The full setup described there is not necessary, but you need to run the commands up to and including
`scripts/setup.sh`.

That is, move to the watching-mb repo and run the following:
```
conda activate bito
pip install -e .
conda env update --file environment.yml

git submodule update --init --recursive
make -C spr_neighbors
conda env config vars set WTCH_ROOT=$PWD
conda activate bito

scripts/setup.sh
```


## Performing an nni-search

For your first try, follow the instructions below and restrict to a top-pruning run on ds3
with uniform branch prior with 50 iterations. That will run in a short amount of time.


### Preparing the data files
Various data files are required to run an nni-search. These files are prepared for the ds-datasets,
assuming you have access to `/fh/fast/matsen_e/shared/vip/ds-golden/`, by running
```
./prep-ds-data.sh
```
This script prepares files based on Mr. Bayes runs with both a uniform or exponential branch prior
on ds1, ds3, ds4, ds5, ds6, ds7, and ds8. You may want to edit the first few lines of the shell script
to restrict to a subset of these runs. Prepping all datasets will take a few hours.

The files for the flu100 dataset are copied from 
[vbpi-torch](https://github.com/zcrabbit/vbpi-torch/blob/main/unrooted/data) and already 
present in `data/flu100`.


### Running the actual search
You can run the nni-search on these datasets and Mr Bayes posteriors with the script
```
./run-nni-search.sh
```
This script requires you have prepared the data files with `prep-ds-data.sh`.
Again, you may want to edit the first few lines of the script to restrict to a subset.
Scripts are available to run a single search on a single dataset, such as `tpds1.sh` and `gpds1.sh`.


### Getting similar statistics for MCMC
Comparison stats for short MCMC runs are generated by running the script
```
./prep-mb-data.sh
```
This script requires you have already prepared the data files with `prep-ds-data.sh`.

### Getting stats for the empirical posterior
Additional stats (used for the credible subsplits plots) are generated by running the 
python script: `python posterior_sdag_stats.py`


### Running a search with multiple starting trees (from the posterior).
To use the first few highest posterior density trees on the ds-datasets, run `multiple_starts_all_ds.sh`.
For the flu100 data set, run `multiple_starts_flu100.sh`.


### Running a search with multiple starting trees (from RAxML).

Install [RAxML](https://github.com/amkozlov/raxml-ng#installation-instructions), if not 
already installed.

Get the trees from RAxML by excuting each of the `run.sh` scripts for the ds-datasets 
and flu100 in the directory `multiple_trees_data/raxml_tests/`, process the trees (for 
all data sets at once) by calling `multiple_trees_data/count_distinct.py`.
and perform a search with one of the `multiple_raxml_starts_ds#.sh`.
That is, run
```
cd multiple_trees_data/raxml_tests/ds1
./run.sh
cd ../ds3/
./run.sh
cd ../ds4/
./run.sh
cd ../ds5/
./run.sh
cd ../ds6/
./run.sh
cd ../ds7/
./run.sh
cd ../ds8/
./run.sh
cd ./fl100/
./run.sh
cd ..
python count_distinct.py
cd ../..
./mutliple_raxml_starts_ds1.sh
./mutliple_raxml_starts_ds3.sh
./mutliple_raxml_starts_ds4.sh
./mutliple_raxml_starts_ds5.sh
./mutliple_raxml_starts_ds6.sh
./mutliple_raxml_starts_ds7.sh
./mutliple_raxml_starts_ds8.sh
./mutliple_raxml_starts_flu100.sh
```

### Run-time experiments.
Scripts for running a search with stats collection off and only computing run-time, use 
a script named something like `time_tpds1.sh` or `time_gpds8.sh`.







