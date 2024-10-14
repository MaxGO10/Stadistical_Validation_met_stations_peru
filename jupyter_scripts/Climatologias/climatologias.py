#Importar librerias

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

#Importar archivos

path = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\Datos.nc' #Ubicacion del archivo
datos = xr.open_dataset(path) #Importar archivo

#Transformando el archivo .nc a matrices dependiendo su variable

datos_pp_values = datos.tp.values #Matriz 3 dimensiones - precipitacion
datos_pp_values = datos_pp_values[:360, :, :]
datos_t_values = datos.t2m.values #Matriz 3 dimensiones - temperatura ambiente
datos_t_values = datos_t_values[:360, :, :]
datos_td_values = datos.d2m.values #Matriz 3 dimensiones - temperatura rocio
datos_td_values = datos_td_values[:360, :, :]
lat = datos.latitude.values #Matriz 1 dimension - latitud
lon = datos.longitude.values #Matriz 1 dimension - longitud
time = datos.time.values #Matriz 1 dimension - tiempo
time = time[:360]

#Transoformando las matrices a los valores tmean, pp, HR

datos_t_C_values = datos_t_values - 273 #Transformando los valores de °K a °C - temperatura
datos_td_C_values = datos_td_values - 273 #Transformando los valores de °K a °C - temperatura rocio
datos_pp_m_values = datos_pp_values * 1000 #Transformando los valores de m a mm - precipitacion
ea = 6.11 * np.exp((7.5 * datos_td_C_values)/(datos_td_C_values + 237.5)) #Hallando la presion vapor agua 
e = 6.11 * np.exp((7.5 * datos_t_C_values)/(datos_t_C_values + 237.5)) #Hallando la presion de saturación
HR = ea/e #Hallando la humedad relativa

#Tranformando a xarray las matrices

datos_t = xr.Dataset(data_vars={'t': (('T', 'Y', 'X'), datos_t_C_values)}, coords={'T': time, 'Y': lat, 'X': lon}) #xarray para la temperatura
datos_pp = xr.Dataset(data_vars={'pp': (('T', 'Y', 'X'), datos_pp_m_values)}, coords={'T': time, 'Y': lat, 'X': lon}) #xarray para la precipitacion
datos_HR = xr.Dataset(data_vars={'HR': (('T', 'Y', 'X'), HR)}, coords={'T': time, 'Y': lat, 'X': lon}) #xarray para la HR

#Hallando la climatologia en los xarray

datos_t_mean = datos_t.groupby('T.month').mean('T') #Climatologia mensual de temperatura
datos_pp_acu = datos_pp.groupby('T.month').sum('T') #Climatologia mensual de precipitacion
datos_HR_mean = datos_HR.groupby('T.month').mean('T') #Climatologia mensual de humedad relativa

# Combinar los Datasets usando xr.merge
clim = xr.merge([datos_t_mean, datos_pp_acu, datos_HR_mean])

#Guardar como archivo nc

output_directory = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Climatologias' #Define el directorio de salida y el nombre del archivo
file = 'climatologia.nc' #Nombre del archivo
output_filepath = f'{output_directory}/{file}' #Combina la ruta del directorio y el nombre del archivo
clim.to_netcdf(output_filepath) # Guarda el xarray en un archivo NetCDF



