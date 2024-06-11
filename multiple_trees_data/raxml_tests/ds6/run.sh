#!/bin/bash

raxml-ng --redo --msa ds6.fasta --outgroup 6 --model JC --tree rand{200} --seed 8675309

