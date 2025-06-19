#!/bin/bash

# Xin Qiao, Jun 19 2025

for file in hmm/*.tar.gz; do
  tar -xvzf "$file" -C hmm
done
