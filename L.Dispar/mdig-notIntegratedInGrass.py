#!/usr/bin/env python
#  This script will run the mdig simulations without using GRASS-GIS.
# Audrey Lustig
# Bio-Protection Research Centre, Lincoln University, NZ
# March 2016


import random
import math
import numpy
import glob, os
import shutil
import multiprocessing
from time import time as now
import numba
import itertools
import time
import sys
import gzip

@numba.jit
def inc_at(distribution, x, y, deltas, survMap):
	# Inner loop of move compiled with Numba
	max_dist = 0.0
	tot_dist = 0.0
	n_dist = 0
	(width, height) = distribution.shape
	for i in range(deltas.shape[0]):
		(dx, dy) = deltas[i]
		newx = x + dx
		newy = y + dy
		if 0 <= newx and newx < width and 0 <= newy and newy < height and not survMap[newx, newy]:
			distribution[newx, newy] += 1
			dist = math.sqrt(dx**2 + dy**2)
			tot_dist += dist
			n_dist += 1
			if dist > max_dist: max_dist = dist
	return (max_dist, tot_dist, n_dist)

###########################################################################
# define a envrionment
class envir:
 # initializing the environment
	def __init__(self,a, lag_time, nbp, dist, r, K, nb_generation, prob_disp,survMap,density,OA,ROA,ROS,mean_dist,max_dist,average_map,extinction):
		(self.w, self.h) = survMap.shape
		assert survMap.shape == survMap.shape

		#### dispersal parameters	
		self.dist=dist # max dispersal distance (cauchy distribution)

		#### population dynamic paramaters
		self.r = r # reproduction rate
		self.K = K # carrying capacity
		self.a = a # allee effect: threshold population size

		#### Introduction parameters
		self.nbp = nbp # number of individual introduced
		self.lag_time = lag_time # time lag between introduction

		self.nb_generation=nb_generation # generation number (max time step)
		self.prob_disp=prob_disp # probability of dispersal
		
		self.survMap=survMap # survival map or ressource distribution map
		# Randomize initial distribution
		self.initMap = numpy.zeros([self.w, self.h], int)
		self.initMap[numpy.random.random_integers(0,self.w-1),numpy.random.random_integers(0,self.h-1)] =1 
		self.distribution = numpy.zeros([self.w, self.h], int)
		#### initialize the population
		self.unsuitable = numpy.count_nonzero(survMap) # number of unsuitable cell in the landscape
		self.occupied=0 # number of cell occupied
		
		self.introduction()

		#### output parameters
		self.density = density # density
		self.OA=OA	# occupied area
		self.ROA=ROA	# rate occupied area
		self.ROS=ROS	# rate of spread
		self.mean_dist=mean_dist # mean dispersal distance
		self.max_dist=max_dist	# max dispersal distance
		self.average_map=average_map # average distribution
		self.extinction=extinction # prob of extinction

	#### Introductions
	# Introduce nbp individuals in the landscape based on self.initMap
	def introduction(self):
		self.distribution[numpy.nonzero(self.initMap)] += self.nbp

	#### Dispersal
	def move(self):
		max_dist = 0.0
		tot_dist = 0.0
		n_dist = 0
		next_distribution = numpy.zeros_like(self.distribution)
		for x in range(self.w): 
			for y in range(self.h):
				n = self.distribution[x, y]
				if not n: continue
				deltas = numpy.zeros([n, 2], int)
				n_far = numpy.random.binomial(n, self.prob_disp) # Binomial distribution to determine number of long dispersal event
				if (n-n_far and n > 0.85*self.K): # if disperse locally
					deltas[n_far:] = numpy.random.random_integers(-1, 1, [n-n_far, 2]) # local dispersal in one of the adgacent cell
				if n_far: # if dispersal over long-distance
					deltas[:n_far] = self.dist * numpy.random.standard_cauchy([n_far,2]) # Cauchy distribution
				(cell_max_dist, cell_tot_dist, cell_n_dist) = inc_at(next_distribution, x, y, deltas, self.survMap)
				if cell_max_dist > max_dist: max_dist = cell_max_dist
				n_dist += cell_n_dist 
				tot_dist += cell_tot_dist
		new_cells = numpy.logical_and(numpy.logical_not(self.distribution), next_distribution).sum()
		self.distribution = next_distribution
		return (new_cells, max_dist, tot_dist, n_dist)

	#### Population dynamic
	def population_growth(self, n):
		if self.a == 0:
			newpop =  n*numpy.exp(self.r *(1 - n/self.K))   # without Allee Effect
		else:
			newpop =  n*numpy.exp(self.r * (1 - n/self.K) * ((n-self.a)/self.K)) # with Allee effect
		return newpop
	
	def reproduce(self):
		newpop = self.population_growth(self.distribution)
		self.distribution = numpy.floor(newpop + numpy.random.uniform(0.0, 1.0, newpop.shape)).astype(int)

	#### Process map
	# calculate for each time step, the average rate of spread, occupied area and others
	def processMap(self, time, new_cells, max_dist, tot_dist, n_dist):
		# area occupied
        
		nb = numpy.sum(self.distribution)
		ncell = numpy.count_nonzero(self.distribution)
		
		# occupied, relative occupied cells
		o = ncell/(self.w*self.h*1.0)
		r = (ncell/((self.w*self.h-self.unsuitable)*1.0))

		# mean, max average distance
		if n_dist > 0:
			mmean = tot_dist / n_dist 
			mmax = max_dist
		else:
			mmean = 0 
			mmax = 0
		
		# summary table for self.repet repetition
		self.density[time]+=nb
		self.OA[time]+=o
		self.ROA[time]+=r
		self.ROS[time]+=new_cells
		self.mean_dist[time]+=mmean
		self.max_dist[time]+=mmax
		
		if time == 30:
			if nb==0:
				self.extinction+=1  
			self.average_map+=self.distribution



	def run(self):
		time=0
		self.processMap(time, numpy.count_nonzero(self.distribution), 0, 0, 0)

		# proceed simulations for a certain number of steps
		while(time<self.nb_generation):
			time+=1
			self.reproduce() # Growth function
			(new_cells, max_dist, tot_dist, n_dist) = self.move() # dispersal and sruvival function
			self.processMap(time, new_cells, max_dist, tot_dist, n_dist) # analyse spread
			if self.lag_time > 0 and time % self.lag_time == 0: # introduction at different time lag
				self.introduction()
		return [self.density, self.OA,self.ROA,self.ROS,self.mean_dist,self.max_dist,self.average_map,self.extinction]


#### Main function for scenarios 
def readMap(filename,fff): # function to read survival and initial distribution
	memMap = []
	locMap = open(filename, 'r')
	
	for i in range(0,6):
		locMap.readline()

	if fff==1:
		for line in locMap:
			line = line.split(' ')
			#line.remove('')
			line[-1] = line[-1][0]
			memMap.append(line)
	else:	
		for line in locMap:
			line = line.split(' ')
			line = line[0:128]
			memMap.append(line)
	return (numpy.array(memMap) == '1').astype(int) # should be bool, but numba is misbehaving


def run(repetition, folder, f, survMap, a, lag_time, nbp, dist, r, K, nb_generation, prob_disp):
	# variable definition	
	nb_init=nbp	
	(w, h) = survMap.shape # Extent of the study
	average_map=numpy.zeros([w, h], int) # Spread map
	density = [0] * (nb_generation+1)	# density
	OA = [0] * (nb_generation+1)	 # Occupied area (Number of cell occupied)
	ROA = [0] * (nb_generation+1)	# Relative occupied area (number of cell occupied / percentage of suitable area in the landscape)
	ROS = [0] * (nb_generation+1)	 # Rate of spread (Number of new cell occupied at each time step)
	mean_dist = [0] * (nb_generation+1) # Mean dispersal distance
	max_dist = [0] * (nb_generation+1) # Max dispersal distance
	extinction = 0 # Number of extinction

	for k in range(nb_init):
		for rep in range(repetition):
			gigi=envir(a, lag_time, nbp, dist, r, K, nb_generation, prob_disp,survMap,density,OA,ROA,ROS,mean_dist,max_dist,average_map,extinction) # Initiate simulation
			[density,OA,ROA,ROS,mean_dist,max_dist,average_map,extinction]=gigi.run() # Run the simulation

	# Output
	extinction_array=[numpy.nan] * (nb_generation+1)	
	extinction_array[0]=extinction
	output=numpy.column_stack((numpy.array(density)/(nb_init*repetition*1.0),numpy.array(OA)/(nb_init*repetition*1.0),numpy.array(ROA)/(nb_init*repetition*1.0),numpy.array(ROS)/(nb_init*repetition*1.0),numpy.array(mean_dist)/(nb_init*repetition*1.0), numpy.array(max_dist)/(nb_init*repetition*1.0),numpy.array(extinction_array)/(nb_init*repetition*1.0)))

	outputStat= outputFolder + '/' + f +  '/'
	if not os.path.exists(outputStat):
		os.makedirs(outputStat)
	#outputMaps= MapsFolder + '/' + f +  '/'
	#if not os.path.exists(outputMaps):
		#os.makedirs(outputMaps)
	name_output_stat= outputFolder + '/' + f +  '/allee.' + str(a) + '.intro.' + str(lag_time) + '.nbp.' + str(nbp) + '.dist.' + str(dist) + '.repro.' + str(r) + '.K.' + str(K) + '.prob_disp.' + str(prob_disp) + '.gz'
	#name_output_maps= MapsFolder + '/' + f +  '/allee.' + str(a) + '.intro.' + str(lag_time) + '.nbp.' + str(nbp) + '.dist.' + str(dist) + '.repro.' + str(r) + '.K.' + str(K) + '.prob_disp.' + str(prob_disp) + '.gz'
	numpy.savetxt(name_output_stat, output, fmt='%f', delimiter=' ')
	#numpy.savetxt(name_output_maps, average_map/(nb_init*repetition*1.0), fmt='%f', delimiter=' ')

def main(repetition):  
	# parallelism
	n_cpus = 7
	pool = multiprocessing.Pool(n_cpus)
	try:
		folder = sys.argv[1] # The user give path to folder 
		f = sys.argv[2]
		locMap1 = survivalMapInputFolder + folder+ '/' + f  
		survMap = readMap(locMap1,1) # read susrvival layer
		# parameters in brackets indicate Allee effect, average time lag between each introductions, number of individual introduced, median distance of long-distance dispersal events, rate of increase, carrying capacity, number of generations, frequency of long-dispersal events 
		parameter_combinations = itertools.product([2], [0], [5], [3,5], [0.815,1.223], [50], [100], [0.05])
		pool.starmap(run, [(repetition, folder, f, survMap)+params for params in parameter_combinations])
	finally:
		pool.close()


# Define input folder

homeDirectory=sys.argv[3]

survivalMapInputFolder = homeDirectory + '/L.Dispar/Data/SurvivalLayer/'
#initialDistributionFolder = 'Data/Ldispar_introductionSites/'

# define output folder
outputFolder = homeDirectory + '/L.Dispar/Results/analysis/'
if not os.path.exists(outputFolder):
	os.makedirs(outputFolder)
MapsFolder = homeDirectory + '/L.Dispar/Results/occupancy.envelopes/'
if not os.path.exists(MapsFolder):
	os.makedirs(MapsFolder)

print('Printing arg1')
print(sys.argv[1])

print('Printing arg2')
print(sys.argv[2])

start = time.clock()
start2= time.time()

main(500) # define here the number of replications

elapsed = (time.clock() - start)
elapsed2= (time.time() - start2)
print("Begining in second from epoch", start2)
print("Ending in second from epoch", time.time())
print("Time difference in CPU time and wall clock unit")
print(elapsed, elapsed2)
