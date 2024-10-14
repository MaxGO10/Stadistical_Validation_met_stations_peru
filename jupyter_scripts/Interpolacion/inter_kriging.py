#Importar librerias
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import regionmask
import os
import xarray as xr

#Importar archivos interpolados
path_spliging = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\interpolated_spliging.nc' #Archivo nc spliging
spli = xr.open_dataset(path_spliging)

#Tranformando el xarray a una lista de valores
pp_monthly_spli = spli.tmax.values
pp_monthly_spli = pp_monthly_spli[:120, :, :]
tmax_monthly_spli = spli.tmin.values
tmax_monthly_spli = tmax_monthly_spli[:120, :, :]
tmin_monthly_spli = spli.pp.values
tmin_monthly_spli = tmin_monthly_spli[:120, :, :]

tmean_monthly_spli = (tmax_monthly_spli+tmin_monthly_spli)/2
ea = 6.11 * np.exp((7.5 * tmin_monthly_spli)/(tmin_monthly_spli + 237.5))
e = 6.11 * np.exp((7.5 * tmean_monthly_spli)/(tmean_monthly_spli + 237.5))
HR = ea/e #Hallando la humedad relativa

#Creando una lista de tiempo
time = pd.date_range(start='2000-01-31', end='2009-12-31', freq='M')
latitud = np.arange(-18, -11.75, 0.25)
longitud = np.arange(-72, -66.75, 0.25)

datos_pp = xr.Dataset(data_vars={'pp': (('T', 'Y', 'X'), pp_monthly_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura
datos_t = xr.Dataset(data_vars={'t': (('T', 'Y', 'X'), tmean_monthly_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura
datos_HR = xr.Dataset(data_vars={'HR': (('T', 'Y', 'X'), tmean_monthly_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura

# Combinar los Datasets usando xr.merge
inter_spli = xr.merge([datos_t, datos_pp, datos_HR])

output_directory = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Interpolacion' #Define el directorio de salida y el nombre del archivo
file = 'interpolate_spli_new.nc' #Nombre del archivo
output_filepath = f'{output_directory}/{file}' #Combina la ruta del directorio y el nombre del archivo
inter_spli.to_netcdf(output_filepath) # Guarda el xarray en un archivo NetCDF