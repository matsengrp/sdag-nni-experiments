#!/bin/bash

raxml-ng --redo --msa ds3.fasta --outgroup 7 --model JC --tree rand{200} --seed 8675309

