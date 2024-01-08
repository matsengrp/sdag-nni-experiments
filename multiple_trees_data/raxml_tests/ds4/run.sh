#!/bin/bash

raxml-ng --redo --msa ds4.fasta --outgroup 8 --model JC --tree rand{200} --seed 8675309

