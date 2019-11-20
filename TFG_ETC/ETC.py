import numpy as np
import math
import astropy.units as u #lo utilizamos para decirle al programa con qué unidades estamos trabajando
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.coordinates import get_moon
from astropy.time import Time
from astropy.utils import iers
from astropy.table import Table
import pysynphot as S
from matplotlib import pyplot as plt
iers.conf.auto_max_age = None
iers.conf.auto_download = False
import glob
import os
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFont
    import tkFileDialog
    import ttk
else:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as tkFileDialog

import ETC_functions    
import parameters

params = parameters_interface()

# USO DE LA FUNCIÓN
# Argumentos de entrada necesarios:    
    #Fecha de observación en formato "YYYY-MM-DD" (El dia de la puesta de sol) y hora de observación
observation = Time(params['day']) + 1*u.d   #Se suma un día para tomar como referencia de nuestros cálculos las 0:00 de la noche de observación
print('Noche: {}'.format(observation))
observation_hour = params['hour'];   

    #Obsercatorio
#location_UAM  = EarthLocation(lat=40.54*u.deg, lon=-3.69*u.deg, height=720*u.m)             #Observatorio de la UAM (Madrid)
#location_CAHA = EarthLocation(lat=37.22*u.deg, lon=-2.55*u.deg, height=2168*u.m)            #Observatorio +de Calar Alto (Almería)
#location_ROMU = EarthLocation(lat=28.76*u.deg, lon=-17.88*u.deg, height=2326*u.m)           #Observatorio de Roque de los Muchachos (Islas Canarias)
#location_MAKE = EarthLocation(lat=19.83*u.deg, lon=-155.47*u.deg, height=4215*u.m)          #Observatorio de Mauna Kea (Hawaii)
#location_PARA = EarthLocation(lat=-24.63*u.deg, lon=-70.41*u.deg, height=2635*u.m)

#Seleccionamos un observatorio
location_observatory = EarthLocation(lat=params['observatory_latitude']*u.deg, lon=params['observatory_longitude']*u.deg, height=params['observatory_altitude']*u.m)                                          
print('Posicion observatorio: {}'.format(location_observatory))

    #Utc offsets (seleccionar zona horaria)
#utc_offset_UAM = 1 #Madrid
#utc_offset_CAHA = 1 #Almeria
#utc_offset_ROMU = 0 #Canarias
#utc_offset_MAKE = -10 #Hawaii
#utc_offset_PARA = -3 #Chile

#Seleccionamos el horario local
utc_offset = params['utc']*u.h
print('UTC = ', utc_offset)


    # Objeto a observar
target = []
name = "Objeto 1"
coordinates_object = SkyCoord(params['ra_object']*u.deg,params['dec_object']*u.deg)
# Hace falta argumento de entrada de magnitud del objeto
mag=params['magn_V']
# Filtro de observación
f =params['filter']

# Espectro del objeto a observar
objeto = params['spectral_template']

# Especificación telescopio         ??????????????????????????????
effective_area = params['Aeff']
gain = params['gain']

# EJECUCIÓN DE LAS FUNCIONES CON LOS ARGUMETNOS DE ENTRADA:
# Corrección de la hora con el etc offset para tener todas las horas en UTC al ejecutar las funciones
observation_this_night = observation_hour*u.h - utc_offset + observation 

# Insertar datos del objeto observado en la clase:
ra = coordinates_object.ra # Ascensión recta
dec = coordinates_object.dec # Declinación
target=(Target(name, ra, dec))

# Abrir parámetros del filtro elegido
dir_path = os.path.dirname(os.path.realpath(__file__))
q_e = os.path.join(dir_path,'data','quantum_efficency', params['quantum_efficiency']+'.txt')
# Abrir espectro del objeto
filename = os.path.join(dir_path,'data','spectral_templates', params['spectral_template']+'.fits')
sp = S.FileSpectrum(filename)
sp.convert('flam')
sp.convert('angstrom')

# 
wl, t, fl, qe, ms  = catalog (f, sp, q_e, dir_path)
# Calcular extinción de la atmósfera
k, extinction_coef = k_ext(1, wl, params['observatory_altitude'], params['observatory_rayleigh'], params['observatory_aerosol'], params['observatory_ozone'])
# Cálculo de la masa de aire
airmass = X_airmass(target.zenith_angle_rad(observation_this_night, location_observatory))
# Creación objeto con condiciones de entrada en la función
conditions = Conditions(k, airmass, effective_area, mag)
# Cálculo del número de electrones por segundo en la CCD2.6
R_obj = electrons_per_second(wl, t, fl, qe, conditions, dir_path)/gain
# Brillo del cielo teniendo en cuenta posición del objeto y Luna
mu_sky, delta_mu_moon, mu_sky_final, warnings_array = target.brillo_cielo(location_observatory, observation_this_night, ms, extinction_coef)
# Cálculo del número de cuentas del cielo (R_sky)
mag0 = mag_zero(gain,f,sp,q_e,conditions, dir_path)
R_sky = 10**(-0.4*(mu_sky-mag0))/gain

# ALGUNOS OUTPUTS
print('El número de cuentas por segundo que contamos en nuestra CCD, usando el filtro, es de: ',
      R_obj)
print('El número de cuentas por segundo del cielo es: ',
      R_sky)

R_tot = R_obj + R_sky
N_bias = params['Nbias']
N_sat = params['Nsat']

t_sat = (N_sat - N_bias)/R_tot
t_sat = t_sat
print('El tiempo de saturación es :', t_sat/3600,'h')

t_values = np.linspace(0, t_sat, 2000)
S_R = signal_noise(R_obj, R_sky, gain, N_bias, 1, t_values)

plt.plot(t_values/3600, S_R)
plt.xlabel('Tiempo de obsrvación (horas)')
plt.ylabel('Señal/Ruido (adimensional)')
plt.show()

print(warnings_array)

    