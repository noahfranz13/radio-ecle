#!/bin/bash

# wipe the old version of the directory
if [[ -d .otter ]]; then
    rm -r .otter
else
    mkdir .otter
fi
    
# update the otter directory in this directory
wget https://github.com/astro-otter/otterdb/raw/main/.otter.tar.gz
tar -xvzf .otter.tar.gz

# convert the ECLE data to the otter format
python3 ecle-data-to-otter.py

# delete the tar file
rm .otter.tar.gz
