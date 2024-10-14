#Importar librerias

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import xarray as xr

#Importar archivo d estaciones

path_est = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\Estaciones.xlsx'
est = pd.read_excel(path_est, "Worksheet", engine='openpyxl')
est = est.iloc[:,[2,3,4,5,6]]

#Importar archivos de datos

pp_monthly = pd.read_csv(r"D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\pp_monthly.csv")
tmin_monthly = pd.read_csv(r"D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\tmin_monthly.csv")
tmax_monthly = pd.read_csv(r"D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Datos\tmax_monthly.csv")

#Ordenando los archivos

#Precipitacion
pp_monthly["fecha"] = pd.to_datetime(pp_monthly['Unnamed: 0'])
pp_monthly = pp_monthly.set_index("fecha")
pp_monthly = pp_monthly.drop(columns=["Unnamed: 0"])
pp_monthly = pp_monthly.iloc[:-12] # Eliminar las últimas 10 filas

#Temperatura minima
tmin_monthly["fecha"] = pd.to_datetime(tmin_monthly['Unnamed: 0'])
tmin_monthly = tmin_monthly.set_index("fecha")
tmin_monthly = tmin_monthly.drop(columns=["Unnamed: 0"])
tmin_monthly = tmin_monthly.iloc[:-12]

#Temperatura maxima
tmax_monthly["fecha"] = pd.to_datetime(tmax_monthly['Unnamed: 0'])
tmax_monthly = tmax_monthly.set_index("fecha")
tmax_monthly = tmax_monthly.drop(columns=["Unnamed: 0"])
tmax_monthly = tmax_monthly.iloc[:-12]

#Filtracion de datos
filtrado_estaciones_pp = est[est["code"].isin(pp_monthly.columns)]
filtrado_estaciones_pp = filtrado_estaciones_pp.set_index("code").loc[pp_monthly.columns].reset_index()

filtrado_estaciones_tmin = est[est["code"].isin(tmin_monthly.columns)]
filtrado_estaciones_tmin = filtrado_estaciones_tmin.set_index("code").loc[tmin_monthly.columns].reset_index()

filtrado_estaciones_tmax = est[est["code"].isin(tmax_monthly.columns)]
filtrado_estaciones_tmax = filtrado_estaciones_tmax.set_index("code").loc[tmax_monthly.columns].reset_index()

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#Interpolacion

#Crearse las grillas con resolcion de 0.25
grid_x, grid_y = np.mgrid[-72:-67:25j, -18:-12:21j]

#Crearse una lista vacia
pp_monthly_spline = []
tmax_monthly_spline = []
tmin_monthly_spline = []

#Bucle para pp
for i in range(pp_monthly.shape[0]):
  x = filtrado_estaciones_pp["lon"].to_list()
  y = filtrado_estaciones_pp["lat"].to_list()
  z = pp_monthly.iloc[i].to_list()

  #Interpolacion spline con el metodo nearest
  grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')
  pp_monthly_spline.append(grid_z)

#Bucle para tmax
for i in range(tmax_monthly.shape[0]):
  x = filtrado_estaciones_tmax["lon"].to_list()
  y = filtrado_estaciones_tmax["lat"].to_list()
  z = tmax_monthly.iloc[i].to_list()

  #Interpolacion spline con el metodo nearest
  grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')
  tmax_monthly_spline.append(grid_z)

#Bucl para tmin
for i in range(tmin_monthly.shape[0]):
  x = filtrado_estaciones_tmin["lon"].to_list()
  y = filtrado_estaciones_tmin["lat"].to_list()
  z = tmin_monthly.iloc[i].to_list()

  #Interpolacion spline con el metodo nearest
  grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')
  tmin_monthly_spline.append(grid_z)

#Transoformando las matrices a los valores tmean, pp, HR

tmean_monthly_spline = []
HR_monthly_spline = []
for i in range(len(tmin_monthly_spline)):
    tmean = (tmax_monthly_spline[i] + tmin_monthly_spline[i])/2
    tmean_monthly_spline.append(tmean)

for i in range(len(tmean_monthly_spline)):
    ea = 6.11 * np.exp((7.5 * tmin_monthly_spline[i])/(tmin_monthly_spline[i] + 237.5))
    e = 6.11 * np.exp((7.5 * tmean_monthly_spline[i])/(tmean_monthly_spline[i] + 237.5))
    HR = ea/e #Hallando la humedad relativa
    HR_monthly_spline.append(HR)

#Creando xarray

#Creando una lista de tiempo
time = pd.date_range(start='2000-01-31', end='2009-12-31', freq='M')
latitud = np.arange(-18, -11.75, 0.25)
longitud = np.arange(-72, -66.75, 0.25)

datos_pp = xr.Dataset(data_vars={'pp': (('T', 'Y', 'X'), pp_monthly_spline)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura
datos_t = xr.Dataset(data_vars={'t': (('T', 'Y', 'X'), tmean_monthly_spline)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura
datos_HR = xr.Dataset(data_vars={'HR': (('T', 'Y', 'X'), tmean_monthly_spline)}, coords={'T': time, 'Y': latitud, 'X': longitud}) #xarray para la temperatura

# Combinar los Datasets usando xr.merge
inter_spline = xr.merge([datos_t, datos_pp, datos_HR])

output_directory = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Interpolacion' #Define el directorio de salida y el nombre del archivo
file = 'interpolate_spline.nc' #Nombre del archivo
output_filepath = f'{output_directory}/{file}' #Combina la ruta del directorio y el nombre del archivo
inter_spline.to_netcdf(output_filepath) # Guarda el xarray en un archivo NetCDF