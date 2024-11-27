#!/bin/bash

# Usar nproc para comprobar el número de núcleos disponibles (para -j)
# Ejecutar en la carpeta split_results/

parallel -j 16 treetime ancestral --aln {} --tree ../../data/global_rescatado_RELATIVE.nwk --outdir ancestral_results/{.} ::: *.fas
