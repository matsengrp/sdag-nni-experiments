#!/bin/bash

raxml-ng --redo --msa ds8.fasta --outgroup 4 --model JC --tree rand{200} --seed 8675309

