#!/bin/bash

raxml-ng --redo --msa ds5.fasta --outgroup 1 --model JC --tree rand{200} --seed 8675309

