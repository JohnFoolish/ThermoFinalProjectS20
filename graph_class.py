# -*- coding: utf-8 -*-
"""
Basic Graphing function class. Creates a simple class to wrap the data we want
to display
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class graph_data:
    
    def __init__(self, model_type, xlabel, ylabel, zlabel):
        self.model_type = model_type
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.fig1 = plt.figure()
        self.ax1 = self.fig1.add_subplot(111, projection='3d')
        self.count = 0
        
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
        