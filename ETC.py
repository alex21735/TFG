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

'''
class View(tk.Frame):
    count = 0
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        b = tk.Button(self, text="OK", command=self.new_window)
        b.pack(side="top")

    def new_window(self):
        self.count += 1
        id = "New window #%s" % self.count
        window = tk.Toplevel(self)
        label = tk.Label(window, text=id)
        label.pack(side="top", fill="both", padx=10, pady=10)'''
        
    
def obs_select():
    
    location = select_obs.get()
    print(location)
    
    page = 1
    

def obs_custom():
    print('lol')

obs_window = Tk()
obs_window.title('Lugar de observación')

lbl_select_obs = Label(obs_window, text='Estos son los observatorios registrados :', anchor=W)
lbl_select_obs.grid(column = 0, row = 0)
select_obs = Combobox(obs_window)
select_obs['values']=('CAHA','UAM')
select_obs.grid(column = 1, row = 0)


info_obs = Label(obs_window, text='Introducir de forma manual la información sobre el observatorio : ', anchor=W)
info_obs.grid(column = 0, row = 1,columnspan = 3)
lat_obs_txt = Label(obs_window, text='Latitud :', anchor=W)
lat_obs_txt.grid(column = 0, row = 2)
lat_obs = Entry(obs_window)
lat_obs.grid(column = 1, row = 2)
lat_obs.insert(0, '0')
lat_obs_u = Label(obs_window, text = 'deg')
lat_obs_u.grid(column = 2, row = 2)
long_obs_txt = Label(obs_window, text='Longitud :', anchor=W)
long_obs_txt.grid(column = 0, row = 3)
long_obs = Entry(obs_window)
long_obs.grid(column = 1, row = 3)
long_obs.insert(0, '0')
long_obs_u = Label(obs_window, text = 'deg')
long_obs_u.grid(column = 2, row = 3)

#Aqui añadir mas cosas sobre el observatorio

btn_select_obs = Button(obs_window, text = 'Quiero este :)', command = obs_select)
btn_select_obs.grid(column = 2, row = 0)

btn_custom_obs = Button(obs_window, text = 'OK', command = obs_custom)
btn_custom_obs.grid(column = 2, row = 20)

obs_window.mainloop()



    



