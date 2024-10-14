#Importar librerias
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import regionmask
import os
import xarray as xr

#Importar archivos interpolados
path_spline = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Interpolacion\interpolate_spline.nc' #Archivo nc spline
path_clim = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Climatologias\climatologia.nc'
spli = xr.open_dataset(path_spline)
clim = xr.open_dataset(path_clim)

# Cambiando la estructura de la climatologia

clim_enero = clim.sel(month=1)
clim_febrero = clim.sel(month=2)
clim_marzo = clim.sel(month=3)
clim_abril = clim.sel(month=4)
clim_mayo = clim.sel(month=5)
clim_junio = clim.sel(month=6)
clim_julio = clim.sel(month=7)
clim_agosto = clim.sel(month=8)
clim_septiembre = clim.sel(month=9)
clim_octubre = clim.sel(month=10)
clim_noviembre = clim.sel(month=11)
clim_diciembre = clim.sel(month=12)

time = pd.date_range(start='1999-01-31', end='2009-01-01', freq='AS')

# Crear los DataArray para cada mes
clim_enero_new = xr.concat([clim.sel(month=1)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_febrero_new = xr.concat([clim.sel(month=2)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_marzo_new = xr.concat([clim.sel(month=3)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_abril_new = xr.concat([clim.sel(month=4)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_mayo_new = xr.concat([clim.sel(month=5)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_junio_new = xr.concat([clim.sel(month=6)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_julio_new = xr.concat([clim.sel(month=7)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_agosto_new = xr.concat([clim.sel(month=8)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_septiembre_new = xr.concat([clim.sel(month=9)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_octubre_new = xr.concat([clim.sel(month=10)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_noviembre_new = xr.concat([clim.sel(month=11)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
clim_diciembre_new = xr.concat([clim.sel(month=12)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})

clim_list = [clim_enero_new, clim_febrero_new, clim_marzo_new, clim_abril_new, clim_mayo_new, clim_junio_new, clim_julio_new, clim_agosto_new, clim_septiembre_new, clim_octubre_new, clim_noviembre_new, clim_diciembre_new]

# Filtrar spli para cada mes del año
enero_spli = spli.sel(T=spli['T'].dt.month == 1)
febrero_spli = spli.sel(T=spli['T'].dt.month == 2)
marzo_spli = spli.sel(T=spli['T'].dt.month == 3)
abril_spli = spli.sel(T=spli['T'].dt.month == 4)
mayo_spli = spli.sel(T=spli['T'].dt.month == 5)
junio_spli = spli.sel(T=spli['T'].dt.month == 6)
julio_spli = spli.sel(T=spli['T'].dt.month == 7)
agosto_spli = spli.sel(T=spli['T'].dt.month == 8)
septiembre_spli = spli.sel(T=spli['T'].dt.month == 9)
octubre_spli = spli.sel(T=spli['T'].dt.month == 10)
noviembre_spli = spli.sel(T=spli['T'].dt.month == 11)
diciembre_spli = spli.sel(T=spli['T'].dt.month == 12)

# Crear una lista con todos los meses
spli_list = [enero_spli, febrero_spli, marzo_spli, abril_spli, mayo_spli, junio_spli, julio_spli, agosto_spli, septiembre_spli, octubre_spli, noviembre_spli, diciembre_spli]

#//////////////////////////////////////////////////

spli_values_t = []
spli_values_pp = []
spli_values_HR = []

# Transformando a valores
for i in spli_list:
    spli_values_t.append(i.t.values)

for i in spli_list:
    spli_values_pp.append(i.pp.values)

for i in spli_list:
    spli_values_HR.append(i.HR.values)


clim_values_t = []
clim_values_pp = []
clim_values_HR = []

# Tranformando a valores
for i in clim_list:
    clim_values_t.append(i.t.values)

for i in clim_list:
    clim_values_pp.append(i.pp.values)

for i in clim_list:
    clim_values_HR.append(i.HR.values)

#////////////////////////////////////////////////////////////////////////////////

# MSE

MSE_tmean_spli = []
MSE_pp_spli = []
MSE_HR_spli = []

# np.nansum es
for i in range(len(spli_values_t)):
    mse = np.nansum((spli_values_t[i] - clim_values_t[i])**2, axis=0) / 10
    MSE_tmean_spli.append(mse)

for i in range(len(spli_values_pp)):
    mse = np.nansum((spli_values_pp[i] - clim_values_pp[i])**2, axis=0) / 10
    MSE_pp_spli.append(mse)

for i in range(len(spli_values_HR)):
    mse = np.nansum((spli_values_HR[i] - clim_values_HR[i])**2, axis=0) / 10
    MSE_HR_spli.append(mse)

#////////////////////////////////////////////////////////////////////////////////

# RMSE

# Convertir la lista a un arreglo de NumPy
MSE_tmean_array_spli = np.array(MSE_tmean_spli)
MSE_pp_array_spli = np.array(MSE_pp_spli)
MSE_HR_array_spli = np.array(MSE_HR_spli)

# Calcular la raíz cuadrada de cada elemento del arreglo
RMSE_tmean_spli = np.sqrt(MSE_tmean_array_spli)
RMSE_pp_spli = np.sqrt(MSE_pp_array_spli)
RMSE_HR_spli = np.sqrt(MSE_HR_array_spli)

# Convertir el arreglo de numpy devuelta a una lista
# RMSE_tmean = RMSE_tmean_array.tolist()
# RMSE_pp = RMSE_pp_array.tolist()
# RMSE_HR = RMSE_HR_array.tolist()


#////////////////////////////////////////////////////////////////////////////////
# MAE
MAE_tmean_spli = []
MAE_pp_spli = []
MAE_HR_spli = []

for i in range(len(spli_values_t)):
    mae = np.abs(np.nansum(spli_values_t[i] - clim_values_t[i], axis=0)) / 10
    MAE_tmean_spli.append(mae)

for i in range(len(spli_values_pp)):
    mae = np.abs(np.nansum(spli_values_pp[i] - clim_values_pp[i], axis=0)) / 10
    MAE_pp_spli.append(mae)

for i in range(len(spli_values_HR)):
    mae = np.abs(np.nansum(spli_values_HR[i] - clim_values_HR[i], axis=0)) / 10
    MAE_HR_spli.append(mae)

#//////////////////////////////////////////////////////////////////////////////////

#MAD

#Hallando la climatologia de los datos
spli_mean = spli.groupby('T.month').mean('T')

#Multiplicando por el numero de años al que se analiza
# Crear los DataArray para cada mes
spli_mean_enero_new = xr.concat([spli_mean.sel(month=1)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_febrero_new = xr.concat([spli_mean.sel(month=2)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_marzo_new = xr.concat([spli_mean.sel(month=3)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_abril_new = xr.concat([spli_mean.sel(month=4)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_mayo_new = xr.concat([spli_mean.sel(month=5)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_junio_new = xr.concat([spli_mean.sel(month=6)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_julio_new = xr.concat([spli_mean.sel(month=7)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_agosto_new = xr.concat([spli_mean.sel(month=8)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_septiembre_new = xr.concat([spli_mean.sel(month=9)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_octubre_new = xr.concat([spli_mean.sel(month=10)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_noviembre_new = xr.concat([spli_mean.sel(month=11)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})
spli_mean_diciembre_new = xr.concat([spli_mean.sel(month=12)]*len(time), dim='month').assign_coords(month=time).rename({'month': 'T'})

spli_mean_list = [spli_mean_enero_new, spli_mean_febrero_new, spli_mean_marzo_new, spli_mean_abril_new, spli_mean_mayo_new, spli_mean_junio_new, spli_mean_julio_new, spli_mean_agosto_new, spli_mean_septiembre_new, spli_mean_octubre_new, spli_mean_noviembre_new, spli_mean_diciembre_new]

spli_mean_values_t = []
spli_mean_values_pp = []
spli_mean_values_HR = []

# Tranformando a valores
for i in spli_mean_list:
    spli_mean_values_t.append(i.t.values)

for i in spli_mean_list:
    spli_mean_values_pp.append(i.pp.values)

for i in spli_mean_list:
    spli_mean_values_HR.append(i.pp.values)


#Parametro del MAD
MAD_tmean_spli =[]
MAD_pp_spli = []
MAD_HR_spli = []

for i in range(len(spli_values_t)):
    mad = np.abs(np.nanmedian(spli_values_t[i] - clim_values_t[i], axis=0))
    MAD_tmean_spli.append(mad)

for i in range(len(spli_values_pp)):
    mad = np.abs(np.nanmedian(spli_values_pp[i] - clim_values_t[i], axis=0))
    MAD_pp_spli.append(mad)

for i in range(len(spli_values_HR)):
    mad = np.abs(np.nanmedian(spli_values_HR[i] - clim_values_t[i], axis=0))
    MAD_HR_spli.append(mad)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Creando una lista de coordenadas
time = np.arange(1, 13)
latitud = np.arange(-18, -11.75, 0.25)
longitud = np.arange(-72, -66.75, 0.25)

# Convirtiendo a netcdf - MSE
datos_MSE_tmean_spli = xr.Dataset(data_vars={'MSE_tmean': (('T', 'Y', 'X'), MSE_tmean_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MSE_pp_spli = xr.Dataset(data_vars={'MSE_pp': (('T', 'Y', 'X'), MSE_pp_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MSE_HR_spli = xr.Dataset(data_vars={'MSE_HR': (('T', 'Y', 'X'), MSE_HR_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})

# Convirtiendo a netcdf - RMSE
datos_RMSE_tmean_spli = xr.Dataset(data_vars={'RMSE_tmean': (('T', 'Y', 'X'), RMSE_tmean_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_RMSE_pp_spli = xr.Dataset(data_vars={'RMSE_pp': (('T', 'Y', 'X'), RMSE_pp_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_RMSE_HR_spli = xr.Dataset(data_vars={'RMSE_HR': (('T', 'Y', 'X'), RMSE_HR_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})

# Convirtiendo a netcdf - MAE
datos_MAE_tmean_spli = xr.Dataset(data_vars={'MAE_tmean': (('T', 'Y', 'X'), MAE_tmean_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MAE_pp_spli = xr.Dataset(data_vars={'MAE_pp': (('T', 'Y', 'X'), MAE_pp_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MAE_HR_spli = xr.Dataset(data_vars={'MAE_HR': (('T', 'Y', 'X'), MAE_HR_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})

# Convirtiendo a netcdf - MAD
datos_MAD_tmean_spli = xr.Dataset(data_vars={'MAD_tmean': (('T', 'Y', 'X'), MAD_tmean_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MAD_pp_spli = xr.Dataset(data_vars={'MAD_pp': (('T', 'Y', 'X'), MAD_pp_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})
datos_MAD_HR_spli = xr.Dataset(data_vars={'MAD_HR': (('T', 'Y', 'X'), MAD_HR_spli)}, coords={'T': time, 'Y': latitud, 'X': longitud})


# Combinar los Datasets usando xr.merge
MSE_xarray_spli = xr.merge([datos_MSE_tmean_spli, datos_MSE_pp_spli, datos_MSE_HR_spli])
RMSE_xarray_spli = xr.merge([datos_RMSE_tmean_spli, datos_RMSE_pp_spli, datos_RMSE_HR_spli])
MAE_xarray_spli = xr.merge([datos_MAE_tmean_spli, datos_MAE_pp_spli, datos_MAE_HR_spli])
MAD_xarray_spli = xr.merge([datos_MAD_tmean_spli, datos_MAD_pp_spli, datos_MAD_HR_spli])

output_directory = r'D:\CICLO 2024-1\Geomatica\Practica\Informe_2_oficial\Estadistica'  # Define el directorio de salida y el nombre del archivo
file_1 = 'MSE_spli.nc'  # Nombre del archivo
file_2 = 'RMSE_spli.nc'
file_3 = 'MAE_spli.nc'
file_4 = 'MAD_spli.nc'

output_filepath_1 = f'{output_directory}/{file_1}'  # Combina la ruta del directorio y el nombre del archivo
output_filepath_2 = f'{output_directory}/{file_2}'  # Combina la ruta del directorio y el nombre del archivo
output_filepath_3 = f'{output_directory}/{file_3}'  # Combina la ruta del directorio y el nombre del archivo
output_filepath_4 = f'{output_directory}/{file_4}'  # Combina la ruta del directorio y el nombre del archivo

MSE_xarray_spli.to_netcdf(output_filepath_1)  # Guarda el xarray en un archivo NetCDF
RMSE_xarray_spli.to_netcdf(output_filepath_2)  # Guarda el xarray en un archivo NetCDF
MAE_xarray_spli.to_netcdf(output_filepath_3)  # Guarda el xarray en un archivo NetCDF
MAD_xarray_spli.to_netcdf(output_filepath_4)  # Guarda el xarray en un archivo NetCDF
