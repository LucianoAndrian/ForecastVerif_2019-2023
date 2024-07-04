#!/bin/bash
ruta="/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/"
ruta_nmme="https://ftp.cpc.ncep.noaa.gov/NMME/prob/netcdf/"

# Fecha actual
current_year=$(date +%Y)
current_month=$(date +%m)

# Encontrar el archivo más reciente
latest_file=$(find $ruta -name 'prate_*.nc' -type f -size +0 -print0 | xargs -0 ls -1t | head -n 1)

echo "$latest_file"

if [ -z "$latest_file" ]; then
    start_year=2018
    start_month=01
else
    # Extraer el año y mes del archivo más reciente
    latest_file_base=$(basename "$latest_file")
    latest_year=$(echo "$latest_file_base" | cut -c 7-10)
    latest_month=$(echo "$latest_file_base" | cut -c 11-12)

    # Calcular el siguiente mes
    start_year=$latest_year
    start_month=$((10#$latest_month + 1))
    if [ $start_month -gt 12 ]; then
        start_month=01
        start_year=$((start_year + 1))
    fi

    # Convertir el mes a formato de dos dígitos
    start_month=$(printf "%02d" $start_month)
fi

# Loop sobre los años y meses para descargar los archivos faltantes
year=$start_year
month=$start_month

echo "$month"
while [ "$year" -lt "$current_year" ] || { [ "$year" -eq "$current_year" ] && [ "$month" -le "$current_month" ]; }; do

    echo test	
    mespr=$(printf "%02d" $month)
    file_pp=${ruta}prate_${year}${mespr}_prob.adj.sea.nc
    file_t=${ruta}tmp2m_${year}${mespr}_prob.adj.sea.nc

    # Descargar archivos solo si no existen
    wget --no-cache -O "$file_pp" "${ruta_nmme}prate.${year}${mespr}.prob.adj.seas.nc"
   
    wget --no-cache -O "$file_t" "${ruta_nmme}tmp2m.${year}${mespr}.prob.adj.seas.nc"
    
    # Incrementar el mes y ajustar el año si es necesario
    month=$((10#$month + 1))
    if [ $month -gt 12 ]; then
        month=01
        year=$((year + 1))
    fi
done

