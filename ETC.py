# -*- coding: utf-8 -*-

import numpy as np
import math
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle,get_sun,get_moon
from astropy.utils import iers
iers.conf.auto_max_age = None
iers.conf.auto_download = False
from matplotlib import pyplot as plt
from matplotlib.widgets import TextBox
from tkinter import *
from tkinter.ttk import *

def obs_select():
    print('a')
    
    

root = Tk()
root.title('Calculadora de tiempo de observación')
inputs = Notebook(root)
inputs.pack()

tab1 = Frame(inputs)
inputs.add(tab1, text='Lugar de observación')

lbl_select_obs = Label(tab1, text='Estos son los observatorios registrados :', anchor=W)
lbl_select_obs.grid(column = 0, row = 0)
select_obs = Combobox(tab1)
select_obs['values']=('CAHA','UAM')
select_obs.grid(column = 1, row = 0)


info_obs = Label(tab1, text='Introducir de forma manual la información sobre el observatorio : ', anchor=W)
info_obs.grid(column = 0, row = 1,columnspan = 3)
lat_obs_txt = Label(tab1, text='Latitud :', anchor=W)
lat_obs_txt.grid(column = 0, row = 2)
lat_obs = Entry(tab1)
lat_obs.grid(column = 1, row = 2)
lat_obs.insert(0, '0')
lat_obs_u = Label(tab1, text = 'deg')
lat_obs_u.grid(column = 2, row = 2)
long_obs_txt = Label(tab1, text='Longitud :', anchor=W)
long_obs_txt.grid(column = 0, row = 3)
long_obs = Entry(tab1)
long_obs.grid(column = 1, row = 3)
long_obs.insert(0, '0')
long_obs_u = Label(tab1, text = 'deg')
long_obs_u.grid(column = 2, row = 3)

#Aqui añadir mas cosas sobre el observatorio

btn_select_obs = Button(tab1, text = 'Quiero este :)', command = obs_select)
btn_select_obs.grid(column = 2, row = 0)

tab2 = Frame(inputs)
inputs.add(tab2, text='Objeto a observar')


obs_window.mainloop()