# -*- coding: utf-8 -*-
"""
Basic Graphing function class. Creates a simple class to wrap the data we want
to display
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class graph_data:
    
    def __init__(self, model_type, xlabel, ylabel, zlabel, xrange, yrange):
        self.model_type = model_type
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.xrange = xrange
        self.yrange = yrange
        self.var1 = ""
        self.var2 = ""
        self.fig1 = plt.figure()
        self.ax1 = self.fig1.add_subplot(111, projection='3d')
        self.count = 0
        
        name = xlabel.upper()
        if (name == "T" or name == "TEMPERATURE"):
            self.xlabel = "Temperature [K]"
            self.var1 = "t"
        elif (name == "V" or name == "VOLUME"):
            self.xlabel = "Volume [L]"
            self.var1 = "v"
        elif (name == "P" or name == "PRESSURE"):
            self.xlabel = "Pressure [atm]"
            self.var1 = "p"
        elif (name == "U" or name == "ENERGY" or name == 'E'):
            self.xlabel = "Energy [J]"
            self.var1 = "u"
            
        name = ylabel.upper()
        if (name == "T" or name == "TEMPERATURE"):
            self.ylabel = "Temperature [K]"
            self.var2 = "t"
        elif (name == "V" or name == "VOLUME"):
            self.ylabel = "Volume [L]"
            self.var2 = "v"
        elif (name == "P" or name == "PRESSURE"):
            self.ylabel = "Pressure [atm]"
            self.var2 = "p"
        elif (name == "U" or name == "ENERGY"or name == 'E'):
            self.ylabel = "Energy [J]"            
            self.var2 = "u"
        
        
    def set_model_type(self, new_model_type):
        self.model_type = new_model_type
                
    def plot_point(self, xpoint, ypoint, zpoint, stable):
        
        if (stable):
            point_color = 'green'
            point_label = 'Stable'
        else:
            point_color = 'red'
            point_label = 'Unstable'
            self.count = self.count - 1
            
        if (self.count <= 0):           
            self.scatter = self.ax1.scatter(xpoint, ypoint, zpoint, c = point_color, label = point_label)
            self.count = self.count + 1
        else:
            self.scatter = self.ax1.scatter(xpoint, ypoint, zpoint, c = point_color)


    def display(self):
        
        self.ax1.set_xlabel(self.xlabel)
        self.ax1.set_ylabel(self.ylabel)
        self.ax1.set_zlabel(self.zlabel)
        
        self.ax1.legend()
        plt.title(self.model_type)

        plt.show()
        