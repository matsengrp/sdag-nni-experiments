#!/bin/bash

raxml-ng --redo --msa flu100.fasta --outgroup 1 --model JC --tree rand{200} --seed 8675309

