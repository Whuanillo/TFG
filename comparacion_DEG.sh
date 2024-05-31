#!/usr/bin/env bash

cd
if [ ! -d data_Hunter ]; then
    echo "El directorio data_Hunter no existe. Creando el directorio..."
    mkdir data_Hunter
    cd data_Hunter
else
    echo "El directorio data_Hunter ya existe en el directorio actual."
    cd data_Hunter
fi

if [ ! -e CountsDE.csv ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/Whuanillo/TFG/raw/main/CountsDE.csv
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

if [ ! -e GSE200621_counts.csv ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/Whuanillo/TFG/raw/main/GSE200621_counts.csv
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

if [ ! -e sample_treat.csv ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/Whuanillo/TFG/raw/main/sample_treat.csv
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

if [ ! -e genes_for_HP_0000510.txt ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/Whuanillo/TFG/raw/main/genes_for_HP_0000510.txt
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

if [ ! -e script_diagramas_venn.sh ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/Whuanillo/TFG/raw/main/script_diagramas_venn.sh
    chmod 755 script_diagramas_venn.sh
    ./script_diagramas_venn.sh
else
    echo "El archivo deseado ya existe en el directorio actual."
    chmod 755 script_diagramas_venn.sh
    ./script_diagramas_venn.sh
fi
