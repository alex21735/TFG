#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:44:28 2019

@author: yago
"""

import numpy as np
import os
import glob
from astropy.table import Table

# Tkinter is for python 2; tkinter is for python 3
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


class Parameter(object):
    
        def __init__(self):
            self.widget = None
    
        def get_value(self):
            return None
    
    
class Numeric_parameter(Parameter):

    def __init__(self, name, description,
                 default=np.NaN, min=np.NaN, max=np.NaN):
        self.name = name
        self.description = description
        self.default = default
        self.value = default
        self.min = min
        self.max = max
        self.widget = None

    def create_tk_widget(self, master, row):
        ttk.Label(master, text=self.name+' :').grid(column=0, row=row)
        self.widget = ttk.Entry(master)
        self.widget.grid(column=1, row=row)
        self.widget.insert(0, str(self.value))
        ttk.Label(master, text=self.description).grid(column=2, row=row)

    def update(self):
        number = float(self.widget.get())
        if (number < self.min) | (number > self.max):
            return False
        else:
            self.value = number
            return True

    def set_value(self, value):
        self.value = float(value)
        if self.widget is not None:
            self.widget.delete(0, tk.END)
            self.widget.insert(0, value)


class String_parameter(Parameter):

    def __init__(self, name, description, default=""):
        self.name = name
        self.description = description
        self.default = default
        self.value = default
        self.widget = None

    def create_tk_widget(self, master, row):
        tk.Label(master, text=self.name+' :').grid(column=0, row=row)
        self.widget = tk.Entry(master)
        self.widget.grid(column=1, row=row)
        self.widget.insert(0, str(self.value))
        tk.Label(master, text=self.description).grid(column=2, row=row)

    def update(self):
        self.value = self.widget.get()
        return True

    def set_value(self, value):
        self.value = str(value)
        if self.widget is not None:
            self.widget.delete(0, tk.END)
            self.widget.insert(0, value)


class Option_parameter(Parameter):

    def __init__(self, name, description, options=[]):
        self.name = name
        self.description = description
        self.options = options
        self.value = None
        self.widget = None

    def create_tk_widget(self, master, row):
        tk.Label(master, text=self.name+' :').grid(column=0, row=row)
        self.widget = ttk.Combobox(master)
        self.widget['values'] = self.options
        self.widget.grid(column=1, row=row)
        tk.Label(master, text=self.description).grid(column=2, row=row)

    def update(self):
        self.value = self.widget.get()
        return True

    def set_value(self, value):  # TO DO: test this code
        if value in self.options:
            self.value = value
            if self.widget is not None:
                self.widget.set(value)

    def get_file(self, files):
        return self.value  # ¿o el elemento de files asociado?


class Input_block(object):

    def __init__(self):
        self.params = {}

    def add_parameter(self, name, param):
        self.params[name] = param

    def read_file(self, filename):
        table = Table.read(filename)
        for row in table:
            self.params[row['param']].set_value(row['value'])

    def create_tk_form(self, master, current_row):
        for param in self.params.values():
            param.create_tk_widget(master, current_row)
            current_row += 1


def parameters_interface ():
    
    block_list = []
    
    # %% Observatory
    
    observatory_dict = {}
    for filename in glob.glob('data/observatories/*csv'):
        observatory_dict[os.path.basename(filename)[:-4]] = filename
    
    observatory = Input_block()
    observatory.add_parameter('observatory_name', String_parameter('Nombre', 'Nombre del observatorio', default='Observatorio'))
    observatory.add_parameter('observatory_latitude', Numeric_parameter('Latitud', 'Latitud en grados', min=-90, max=90,default=0))
    observatory.add_parameter('observatory_longitude', Numeric_parameter('Longitud', 'Longitud en grados', min=-180, max=180,default=0))
    observatory.add_parameter('observatory_altitude', Numeric_parameter('Altitud', 'Altitud sobre el nivel del mar en metros', min=0,default=0))
    observatory.add_parameter('observatory_rayleigh', Numeric_parameter('Exitinción Rayleigh', 'coeficiente de dispersión para el observatorio', min=0,default=1))
    observatory.add_parameter('observatory_aerosol', Numeric_parameter('Extinción por aerosoles', 'coeficiente de dispersión por aerosoles para el observatorio', min=0,default=0))
    observatory.add_parameter('observatory_ozone', Numeric_parameter('Extinción por ozono', 'coeficiente de absorción del ozono para el observatorio', min=0,default=0))
    
    block_list.append(observatory)
    
    
    # %% Instrument
    
    instrument_dict = {}
    for filename in glob.glob('data/instruments/*csv'):
        instrument_dict[os.path.basename(filename)[:-4]] = filename
    
    filter_dict = {}
    for filename in glob.glob('data/filters/*txt'):
        filter_dict[os.path.basename(filename)[:-4]] = filename
        
    qe_dict = {}
    for filename in glob.glob('data/quantum_efficency/*txt'):
        qe_dict[os.path.basename(filename)[:-4]] = filename
    
    instrument = Input_block()
    instrument.add_parameter('instrument_name', String_parameter('Nombre', 'Nombre del instrumento',default='CCD'))
    instrument.add_parameter('Aeff', Numeric_parameter('Area efectiva', 'Área efectiva, en m^2', min = 0,default=1))
    instrument.add_parameter('gain', Numeric_parameter('Gain', 'Ganancia inversa, en e^-/ADU', min = 0,default=1))
    instrument.add_parameter('pixel_size', Numeric_parameter('Tamaño del pixel', 'Tamaño del pixel en arcsec', min = 0,default=1))
    instrument.add_parameter('Nbias', Numeric_parameter('N bias', 'Numero de cuentas de bias por segundo', min = 0,default=300))
    instrument.add_parameter('Nsat', Numeric_parameter('N saturacion', 'Numero de cuentas máximo que podemos almacenar', min = 0,default=65536))
    instrument.add_parameter('quantum_efficiency',
                             Option_parameter('Eficiencia cuántica', 'en función de la longitud de onda',
                                              options=list(qe_dict.keys())))
    instrument.add_parameter('filter',
                             Option_parameter('Filtro', 'Filtro fotométrico',
                                              options=list(filter_dict.keys())))
    
    block_list.append(instrument)
    
    
    # %% Source
    
    spectral_template_dict = {}
    for filename in glob.glob('data/spectral_templates/*fits'):
        spectral_template_dict[os.path.basename(filename)[:-5]] = filename
    
    source = Input_block()
    source.add_parameter('ra_object', Numeric_parameter('RA', 'Ascensión Recta en grados', min=0, max=360, default=0))
    source.add_parameter('dec_object', Numeric_parameter('DEC', 'Declinación en grados', min=-90, max=90, default=0))
    source.add_parameter('magn_V', Numeric_parameter('V', 'Magnitud aparente en banda V en el sistema de referencia de Vega', default=22.4))
    source.add_parameter('spectral_template',
                         Option_parameter('Espectro', 'Forma aproximada de la distribución espectral de energía',
                                          options=list(spectral_template_dict.keys())))
    
    block_list.append(source)
    
    
    # %% Observing conditions
    
    observing_conditions = Input_block()
    observing_conditions.add_parameter('day', String_parameter('Día de observación', ' en formato YYYY-MM-DD ',default='2020-01-01'))
    observing_conditions.add_parameter('hour', Numeric_parameter('Hora de observación', 'hora local', min=0, max=23, default=0))
    observing_conditions.add_parameter('utc', Numeric_parameter('Huso horario', ' en el lugar de observación', min=-14, max=14, default=0))
    
    block_list.append(observing_conditions)
    
    
    # %% Create GUI window
    
    root = tk.Tk()
    root.title('Calculadora de tiempo de observación')
    inputs = ttk.Notebook(root)
    
    # Observatory
    
    tab1 = ttk.Frame(inputs)
    inputs.add(tab1, text='Observatorio')
    observatory_selector = ttk.Combobox(tab1)
    observatory_selector.grid(row=0)
    observatory_selector['values'] = list(observatory_dict.keys())
    observatory_selector.bind(
            '<<ComboboxSelected>>',
            lambda x: observatory.read_file(observatory_dict[
                    observatory_selector.get()]))
    observatory.create_tk_form(tab1, 1)
    
    # Instrument
    
    tab2 = ttk.Frame(inputs)
    inputs.add(tab2, text='Instrumento')
    instrument_selector = ttk.Combobox(tab2)
    instrument_selector.grid(row=0)
    instrument_selector['values'] = list(instrument_dict.keys())
    instrument_selector.bind(
            '<<ComboboxSelected>>',
            lambda x: instrument.read_file(instrument_dict[
                    instrument_selector.get()]))
    instrument.create_tk_form(tab2, 1)
    
    # Source
    
    tab3 = ttk.Frame(inputs)
    inputs.add(tab3, text='Fuente')
    source.create_tk_form(tab3, 1)
    
    # Observing conditions
    
    tab4 = ttk.Frame(inputs)
    inputs.add(tab4, text='Observación')
    observing_conditions.create_tk_form(tab4, 1)
    
    # Footer
    
    inputs.pack()
    inputs.grid(row=0)
    
    
    parameters = {}
    def validate():  # TO DO: Validation ;^D
        for block in block_list:
            for param in block.params:
                block.params[param].update()
                parameters[param] = block.params[param].value
        root.destroy()
    
    
    tk.Button(root, text='OK', command=validate).grid(row=1)
    root.mainloop()
    
    
    # %% Return parameters dictionary
    
    return parameters

# -----------------------------------------------------------------------------
#                                                    ... Paranoy@ Rulz! ;^D
# -----------------------------------------------------------------------------
