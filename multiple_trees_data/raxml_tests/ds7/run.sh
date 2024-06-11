#!/bin/bash

raxml-ng --redo --msa ds7.fasta --outgroup 7 --model JC --tree rand{200} --seed 8675309

