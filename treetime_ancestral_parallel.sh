#!/bin/bash

# Usar nproc para comprobar el número de núcleos disponibles (para -j)
# Ejecutar en la carpeta split_results/

parallel -j 3 treetime ancestral --aln {} --tree ../../data/global_rescatado_RELATIVE_CLEANED.nwk --outdir ancestral_results/{.} ::: *.fas
