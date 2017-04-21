#!/bin/bash

#SBATCH -J checkJob

file="$1"

if [ -f $file ]; then
	rm $file
fi

echo "COMPLETED $1" > $file
