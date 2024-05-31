#!/usr/bin bash
pwd
# GENERACIÓN DEL TRANSCRIPTOMA DE REFERENCIA - OPTIMIZADO

#Creamos el directorio donde vamos a crear el TRANSCRIPTOMA en caso de que no exista
cd
if [ ! -d db ]; then
    echo "El directorio db no existe. Creando el directorio..."
    mkdir db
    cd db
    mkdir refanno
    cd refanno
else
    echo "El directorio db ya existe en el directorio actual."
    cd ~/db/refanno
fi

#Obtenemos la secuencia del genoma
if [ ! -e GRCm39.primary_assembly.genome.fa ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/GRCm39.primary_assembly.genome.fa.gz
    gunzip GRCm39.primary_assembly.genome.fa.gz
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Obtenermos las secuencias de transcripciones
if [ ! -e gencode.vM34.transcripts.fa ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.transcripts.fa.gz
    gunzip gencode.vM34.transcripts.fa.gz
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Obtenemos anotaciones del genoma
if [ ! -e *.gtf ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz
    gunzip gencode.vM34.annotation.gtf.gz
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Limpiamos el cabezero del fasta (fasta header)
if [ ! -e *.clean.fa ]; then
    echo "El archivo deseado no existe. Descargando..."
    sed 's/|ENSG.*//' gencode.vM34.transcripts.fa > gencode.vM34.transcripts.clean.fa
    grep ">" gencode.vM34.transcripts.fa | head
    grep ">" gencode.vM34.transcripts.clean.fa | head
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Creamos informacion tabular desde GTF
if [ ! -e *.txt ]; then
    echo "El archivo deseado no existe. Descargando..."
    cat gencode.vM34.annotation.gtf | \
    awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' | \
    sed 's/gene_id "//' | \
    sed 's/gene_name "//'| \
    sed 's/gene_type "//' | \
    sed 's/"//g' | \
    awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
    sort -k3,3 -k4,4n -k5,5n | \
    sed "1i\GeneID\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
    > gencode.vM34.annotation_genes.txt

    cat gencode.vM34.annotation.gtf | \
    awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"transcript") print a[1]"\t"a[4]"\t"$1":"$4"-"$5"\t"a[2]"\t"a[6]"\t"a[5]"\t"$7}' | \
    sed 's/gene_id "//' | \
    sed 's/gene_name "//'| \
    sed 's/transcript_id "//' | \
    sed 's/transcript_name "//' | \
    sed 's/transcript_type "//' | \
    sed 's/"//g' | \
    awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"$6"\t"$7"\t"a[3]-a[2];}' | \
    sort -k3,3 -k4,4n -k5,5n | \
    sed "1i\GeneID\tGeneSymbol\tChromosome\tStart\tEnd\tTranscriptID\tTranscriptName\tClass\tStrand\tLength" \
    > gencode.vM34.annotation_transcripts.txt
else
    echo "Los archivos deseados ya existen en el directorio actual."
fi

#Tras esto descargamos herramientas que utilizaremos
cd
if [ ! -d tools ]; then
    echo "El directorio db no existe. Creando el directorio..."
    mkdir tools
    cd tools
else
    echo "El directorio db ya existe en el directorio actual."
    cd tools
fi

#Descarga de FastQC
if [ ! -e fastqc_v0.11.9.zip ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    unzip fastqc_v0.11.9.zip
    cd FastQC
    chmod 755 fastqc
    ./fastqc --help
    cd ..
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Descarga de BBTools
if [ ! -e BBMap_38.82.tar.gz ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget -O BBMap_38.82.tar.gz 'https://downloads.sourceforge.net/project/bbmap/BBMap_38.82.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2FBBMap_38.82.tar.gz%2Fdownload&ts=1587572478'
    tar xvf BBMap_38.82.tar.gz
    cd bbmap
    ./bbduk.sh
    cd ..
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Descarga de Salmon
if [ ! -e salmon-1.2.1_linux_x86_64.tar.gz ]; then
    echo "El archivo deseado no existe. Descargando..."
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.2.1/salmon-1.2.1_linux_x86_64.tar.gz
    tar xvf salmon-1.2.1_linux_x86_64.tar.gz
    mv salmon-latest_linux_x86_64 salmon-1.2.1
    cd salmon-1.2.1/bin
    ./salmon 
    cd ~/tools
else
    echo "El archivo deseado ya existe en el directorio actual."
fi

#Preparamos los indexados
cd ~/db/refanno

if [ ! -d gencode.vM34_salmon-1.2.1 ]; then
    echo "El directorio deseado no existe. Creando el directorio..."
    ~/tools/salmon-1.2.1/bin/salmon index -p 6 --gencode -t gencode.vM34.transcripts.fa -i gencode.vM34_salmon-1.2.1
else
    echo "El directorio ya existe en el directorio actual."
fi
