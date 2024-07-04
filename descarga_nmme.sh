#!bin/bash
#download nmme operational forecast
ruta="/pikachu/datos/luciano.andrian/verif_2019_2023/nmme_pronos/"
ruta_nmme="https://ftp.cpc.ncep.noaa.gov/NMME/prob/netcdf/"

for k in {2024..2024}; do
	for l in `seq -f "%02g" 1 12`; do
		FM=`expr $[10#$l]`
		FY=$k
		
		if [ $((10#$FM)) -le 9 ] ; then
			mespr=0${FM#0}
		else
			mespr=$FM
		fi

		file_pp=${ruta}prate_${k}${mespr}_prob.adj.sea.nc
		file_t=${ruta}tmp2m_${k}${mespr}_prob.adj.sea.nc

		wget --no-cache -O "$file_pp" "${ruta_nmme}prate.${k}${mespr}.prob.adj.seas.nc"
		wget --no-cache -O "$file_t" "${ruta_nmme}tmp2m.${k}${mespr}.prob.adj.seas.nc"
done
done


