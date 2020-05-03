import sys
import math
import numpy as np
from matplotlib import pyplot as plt

"""python code reading log.lammps file of nvt_BB_real, create 
	tempeture graph as a function of time"""

	
	
"""create tempeture graph using mathplot
	input- axis x, axis y, label for each axis and title for the graph"""
def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph", name="graph.png"):
	plt.plot(axis_x,axis_y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.savefig(name)
	plt.show()
	

"""create Graphs of tempeture function of time"""		
def log_graphs():
	fp = open("log.lammps_nvt.txt","r") 
	text=fp.read() 
	text_list = text.split("\n")
	timeStep_arr = []
	temp_arr = []
	for line in text_list[66:-1]:
		timeStep_arr.append(int(line.split()[0]))
		temp_arr.append(float(line.split()[1]))

	make_graph(timeStep_arr, temp_arr, "TimeStep", "Temprature", "Temprature As A Function Of Time","temp.png")

		
#operate functions as you want to create the graphs you want.	
log_graphs()