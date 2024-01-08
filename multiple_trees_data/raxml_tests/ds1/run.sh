#!/bin/bash

raxml-ng --redo --msa ds1.fasta --outgroup 15 --model JC --tree rand{200} --seed 8675309

