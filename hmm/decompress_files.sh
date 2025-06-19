#!/bin/bash

# Xin Qiao, Jun 19 2025

for file in ./*.tar.gz; do
  tar -xvzf "$file"
done
