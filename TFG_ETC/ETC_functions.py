#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import pysynphot as S

import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_sun
from astropy.coordinates import get_moon
from astropy.utils import iers
iers.conf.auto_max_age = None
iers.conf.auto_download = False


#  Objeto a observar y funciones dependientes del objeto
class Target:

    def __init__(self, name, RA, DEC):
        # Entrada: nombre de la estrella, RA: ascensión recta (hh:mm:ss), DEC: declinación (dd:mm:ss)
        self.coord = SkyCoord(ra=RA, dec=DEC)
        self.name = name

    def lunar_distance_rad(self, time, moon):
        # Salida: separación angular (rad) entre la luna y el objeto
        moon_data = moon.separation(self.coord)
        return moon_data.rad

    # Ángulo cenital en radianes
    def zenith_angle_rad(self, time, earth_location):
        # Entrada: fecha de observación, localización del observatorio, objeto
        # Salida: ángulo cenital en grados, minutos y segundos
        a = self.coord.transform_to(AltAz(obstime=time, location=earth_location))
        z = np.pi/2-a.alt.rad
        return z

    def brillo_cielo(self, location_observatory, time, ms, extinction_coef):
        # Entrada: posición del observatiorio, fecha y hora (UTC), brillo del cielo para
        # el filtro elegido, coeficiente de extinción del cielo (medio para el filtro)
        def lunar_phase(time):
            # Cálculo del ángulo de la fase lunar
            # 0 para luna llena
            # 180 para luna nueva
            sun = get_sun(time)
            moon = get_moon(time)
            elongation = sun.separation(moon)
            lunar_angle = np.arctan2(sun.distance*np.sin(elongation),
                                     moon.distance - sun.distance*np.cos(elongation))
            angle = lunar_angle/u.rad*180/np.pi  # ángulo en grados
            return angle

        # BRILLO DEL CIELO SIN LA LUNA
        # Variación con el ángulo cenital
        zenith = self.zenith_angle_rad(time, location_observatory)
        mz = ms - 2.5*np.log10(X_airmass(zenith)) + extinction_coef*(X_airmass(zenith)-1)

        # Posición de la Luna y el objeto
        iers.conf.auto_download = False
        moon_data = get_moon(time)        # Datos de la Luna
        lunar_angle = lunar_phase(time)   # Ángulo de fase de la luna en grados
        lunar_separation = self.lunar_distance_rad(time, moon_data)      # Distancia angular entre la Luna y el objeto
        a = moon_data.transform_to(AltAz(obstime=time, location=location_observatory))
        moon_zenith_angle_rad = np.pi/2 - a.alt.rad  # Ángulo cenital de la luna
        target_zenith_angle_rad = self.zenith_angle_rad(time, location_observatory)   # Ángulo cenital del objeto observado

        def B_moon_nanoLamberts(lunar_angle, lunar_separation, moon_zenith_angle, target_zenith_angle):
            illuminance = 10**(-0.4*(3.84 + 0.026*np.abs(lunar_angle) + 4e-9*lunar_angle**4))  # Illuminance of the moon inside the atmosphere
            f_rho = 10**(5.36)*(1.06 + (np.cos(lunar_separation))**2) + 10**(6.15-lunar_separation/40)  # Función de scattering (rho es la distancia entre la luna y el objeto)
            B_moon = f_rho * illuminance * 10**(-0.4*extinction_coef
                                                * X_airmass(moon_zenith_angle)
                                                ) * (1 - 10**(-0.4*extinction_coef*X_airmass(target_zenith_angle)))
            return B_moon  # Brillo de la luna en nano lamberts

        B_moon = B_moon_nanoLamberts(lunar_angle, lunar_separation, moon_zenith_angle_rad, target_zenith_angle_rad)
        # Pasamos a nano lamberts el valor obtenido para el cielo sin Luna en magnitudes
        # con la transformación de unidades dada por KRISCIUNAS
        B_sky = 34.08*np.exp(20.7233 - 0.92104 * mz)  # Sky brighntess in nano Lamberts
        # Después, calculamos la variación en magnitud que produce el contraste con la luna
        # Al calcular la variación relativa podemos ignorar las unidades y no hacer el cambio inverso de los nano lamberts
        delta_m_moon = -2.5*np.log10((B_sky+B_moon)/B_sky)

        # Cuadro de warnings
        warnings_array = []
        # Luna llena
        if 0 < lunar_angle < 10:
            warnings_array.append("La Luna está prácticamente llena")
            # Objeto muy cerca de la luna
        if -0.175 < lunar_separation < 0.175:
            warnings_array.append("El objeto está muy cerca de la luna")
        # Masa de aire muy alta
        if X_airmass(target_zenith_angle_rad) > 2:
            warnings_array.append("El valor de la masa de aire es muy alto")

        # mz - brillo del cielo tras corregir con el ángulo cenital (magnitudes)
        # delta_m_moon - Corrección del brillo de la luna (magnitudes)
        # mz+delta_m_moon - Suma de ambas contribuciones
        # warnings_array - Qué está saliendo mal al medir
        return mz, delta_m_moon, mz+delta_m_moon, warnings_array


# Para meter en bloque las condiciones en la función que calcula R_obj
class Conditions:

    def __init__(self, extinction, air_mass, effective_area, magnitud):
        self.extinction = extinction
        self.air_mass = air_mass
        self.effective_area = effective_area
        self.magnitud = magnitud


def catalog(s_filter, spectrum, quantum_efficiency, dir_path):
    #Entradas: string con el nombre del filtro
    #          archivo FITS del espectro deseado
    #          string del path del archivo txt con los parámetros de la ccd
    file = os.path.join(dir_path,'data','filters',s_filter+'.txt')
    if  s_filter == 'U':
        magnitud_cielo = 21.98
    elif s_filter == 'V':
        magnitud_cielo = 22.6
    elif s_filter == 'B':
        magnitud_cielo = 21.5
    elif s_filter == 'R':
        magnitud_cielo = 20.6
    elif s_filter == 'I':
        magnitud_cielo = 18.7

    else:
        return print('Error, no has seleccionado un filtro concreto')

    fwlt = np.loadtxt(file)
    fwl = fwlt[:, 0]  # array con las longitudes de onda del filtro
    ft = fwlt[:, 1]   # array con el factor de transmisión del filtro
    swl = spectrum.wave    # array con las longitudes de onda del objeto
    sflux = spectrum.flux  # array el flujo del objeto por longitud de onda
    ccdwlqe = np.loadtxt(quantum_efficiency)
    ccdwl = ccdwlqe[:, 0]  # array con las longitudes de onda de la ccd
    ccdqe = ccdwlqe[:, 1]  # array con la quantum efficiency
    num = len(swl)+len(fwl)+len(ccdwl)
    wlinterp = np.linspace(fwl[0], fwl[-1], num)  # array longitudes de onda interpoladas (Angstrom (\AA))
    tinterp = np.interp(wlinterp, fwl, ft)        # array transmisión interpolada (adimensional)
    fluxinterp = np.interp(wlinterp, swl, sflux)  # array flujo por longitud de onda interpolado(erg s^-1 cm^−2 \AA^−1)
    efficiencyinterp = np.interp(wlinterp, ccdwl, ccdqe)  # array quantum efficiency (tanto por 1)
    # Salidas: longitud de onda, transmisión, flujo y eficiencia interpolada
    return wlinterp, tinterp, fluxinterp, efficiencyinterp, magnitud_cielo


def mag_zero(Gain, fU, sp, q_e, conditions, dir_path):

    wl, t, Fl, qe, ms = catalog(fU, sp, q_e, dir_path)
    # wl=wl*10**-10  #angstroms to meters
    # Fl = Fl*10**7 #Flams to SI units
    # Bisection Method
    magnitud_1 = 10
    magnitud_2 = 50
    tolerance = 1e-10
    contador = 0
    magnitud_new = 20
    conditions.magnitud = magnitud_new
    while (np.abs(electrons_per_second(wl, t, Fl, qe, conditions, dir_path)/Gain - 1)  >= tolerance
           and contador < 100):

        magnitud_new = 0.5*(magnitud_1 + magnitud_2)
        conditions.magnitud = magnitud_new
        value_R = electrons_per_second(wl, t, Fl, qe, conditions, dir_path)/Gain
        if value_R > 1:
            magnitud_1 = 0.5*(magnitud_1 + magnitud_2)
        else:
            magnitud_2 = 0.5*(magnitud_1 + magnitud_2)
        contador += 1
    return magnitud_new



# Masa de aire (en función del ángulo cenital)
def X_airmass(zenith_angle):
    X = (1-0.972*np.sin(zenith_angle)**2)**(-0.5)
    # Importante meter ángulo en radianes
    return X

# COEFICIENTE DE EXTINCIÓN DEL CIELO
def k_ext(secz, wavelength, height, a,b,c): #Función que calcula el coeficiente de extinción en función de la elevación
    #Se le mete como argumento la masa de aire, el intervalo de longitudes de onda, el numero de puntos y el observatorio
    #Alturas de distintos observatorios:
    h_obs = height/1000
    x = wavelength #* u.Angstrom, vector de longitudes de onda

    #Funciones que luego sirven para calcular la extinción:
    B = 0.0095*(np.exp(-(h_obs/7.996))) #teniendo en cuenta que la altura va en km
    b = 0.0414*(np.exp(-(h_obs/1.500)))


   #Contribuciones a la extinción:
    ext_r = B*((x/10000)**(-4)) #extincion por dispersion rayleigh
    ext_aer = b*((x/10000)**(-0.8)) #extincion por aerosoles (contaminación)
    ext_o = 0.4*np.exp(-((x-6000)/1200)**2) #extincion por la capa de ozono

    #Contribución de agua y oxígeno:
    def gaussiana(y0, A, xc, w, x):
        return y0 + A*np.exp(-(x-xc)**2/(2*(w**2)))

    #CONTRIBUCIÓN DEL OXÍGENO------------------------

    #Desde 6833 a 6933 gaussiana oxígeno:
    A = 0.178 #(amplitud)
    xc = 6870 #(centro)
    w = 7 + (25.6)*(secz - 1)

    #Desde 7600 a 7680 gaussiana oxígeno:
    A2 = 0.765
    xc2 = 7605
    w2 = 28 + (25.6)*(secz - 1)

    #CONTRIBUCIÓN DEL AGUA---------------------------

    #Desde 7133 a 7333 gaussiana:
    A1 = 0.07 #(amplitud)
    xc1 = 7267 #(centro)
    w1 = 60 + 3.59*(secz - 1) #(ancho)

    #Offset promedio para agua y oígeno:
    y0 = 0.0

    f = [a,b,c]

    #Suma total ponderada
    k = ((f[0]*ext_r)+(f[2]*ext_o)+(f[1]*ext_aer) + gaussiana(y0, A, xc, w, x) + gaussiana(y0, A1, xc1, w1, x) + gaussiana(y0, A2, xc2, w2, x))*secz

    #Graficamos:
    #plt.plot(x,k)
    #plt.title('Extinción Atmosférica ' + f[3])
    #plt.xlabel('lambda [Angstrom]')
    #plt.ylabel('K [mag]')
    #plt.grid()
    #plt.show()

    # Sacar el coeficiente de extinción medio para cada filtro
    num = 10000 # Número de pasos de la integral
    step = (wavelength[-1]-wavelength[0])/num
    I1= 2*sum(k)
    I1 = step*0.5*(I1 - k[0]-k[-1])
    I2 = 2*sum(np.ones(len(wavelength)))
    I2 = step*0.5*(I2-2)
    k_num = I1/I2

    return k, k_num



def electrons_per_second (wavelength, transmission, flux, quantum_efficiency, conditions, dir_path):
    #Los inputs son el filtro a usar ('U','B','V','I'), el espectro de referencia a utilizar
    #El tamaño aparente del objeto (en estereorradianes)
    file = os.path.join(dir_path,'data','filters','V.txt')

    extinction = conditions.extinction
    effective_area = conditions.effective_area
    air_mass = conditions.air_mass
    magnitud = conditions.magnitud

    #iniciamos las constantes que se usan para el cuerpo negro (en uds del SI)
    c  = 3e8 #speed light
    h  = 6.626e-34 #Planck Constant
    kb = 1.38e-23

    wavelength = wavelength *(10**-10) #transform angstroms into meters
    #plt.plot(wavelength, transmission)

    #plt.xlabel('Longitud de onda (en m)')
    #plt.ylabel('Transmision (adimensional)')
    #plt.show()

    vega_flux = S.Vega.flux *(10**7)
    vega_wl = S.Vega.wave *(10**-10)
    F0 = np.interp(wavelength, vega_wl, vega_flux)
    #Vamos a calcular ahora las integrales
    a = wavelength[0]
    b = wavelength[-1]
    step = (b-a)/len(wavelength)

    t0file = np.loadtxt(file)
    fwl0 = t0file[:,0] #array con longitudes de onda de V
    ft0 = t0file[:,1]  #array con el factor de transmisión de V
    transmission0 = np.interp(wavelength, fwl0, ft0) #array transmisión interpolada (adimensional)

    #definimos las funciones a integrar
    f1 = flux * transmission0 * (wavelength/(h*c))
    f2 = F0 * transmission0 * (wavelength/(h*c))

    I1 = 2*sum(f1)                       #método del trapecio
    I1= step*0.5*(I1 - f1[0]-f1[-1])
    #print('La integral 1 vale: ',I1)

    I2 = 2*sum(f2)
    I2= step*0.5*(I2 - f2[0]-f2[-1])
    #print('La integral 2 vale: ',I2)
    #La magnitud, por tanto, viene dada por lo siguiente
    mag = -2.5 * np.log10(I1/I2)
    #print('La magnitud que hemos calculado vale : ',mag)
    #Hallamos la constante de normalización
    Normalization = (10**(-0.4*(magnitud-mag)))
    #Fv = Fv*N
    #print('N vale : ',Normalization)

    #Definimos la integral a resolver
    F = ((wavelength/(c*h)) * Normalization *flux* (10**(-0.4*extinction*air_mass))
    * effective_area * quantum_efficiency*transmission)


    #plt.loglog(wavelength, F)
    #plt.show()
    #plt.plot(wavelength, transmission, 'r:')
    #plt.show()

    R = 2*sum(F)
    R = step*0.5*(R - F[0] - F[-1])

    return R

def signal_noise(R_obj, R_sky, gain, N_bias, random_error, t):
    return((R_obj*t)/np.sqrt(gain*(random_error**2+t*(R_obj+R_sky)+N_bias)))

