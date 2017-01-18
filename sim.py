import matplotlib.pyplot as plt
import numpy as np
import shapes
import time
from matplotlib.animation import FuncAnimation

MAX_PLOT_PATHS = 5

def initElectrons(numElectrons, width, height, vtherm, polygons) :
	tempPosition = []
	tempVelocity = []
	paths = []

	for i in range(numElectrons) :
		# Generate random position
		valid = False
		while not valid :
			attemptPos = np.array([np.random.uniform(0, width),
						np.random.uniform(0, height)])
			valid = not isPointInPolys(polygons, attemptPos)
			if valid :
				tempPosition.append([[attemptPos[0],
							attemptPos[1]]])
		

		angle = np.random.uniform(0,2*np.pi)
		tempVelocity.append([[vtherm*np.cos(angle),
				vtherm*np.sin(angle)]])


	for i in range(min(numElectrons,MAX_PLOT_PATHS)) :
		paths.append([[],[]])

	return np.array(tempPosition), np.array(tempVelocity), paths

# This function checks if a given point lies in any of the given polygons
def isPointInPolys(polygons, point) :
	for polygon in polygons :
		inPoly = True
		for wall in polygon :
			originToPoint = np.array([point[0] - wall[0],point[1] - wall[1]])
			wallVec = np.array([wall[3],-1*wall[2]])
			if np.dot(originToPoint, wallVec) > 0 :
				inPoly = False
				break
		if inPoly : return True
	return False
			
	

def iterate(plt, width, timestep, position, velocity, wallPos, wallVec, wallNorm, path) :

	# Collisions
	w = np.tile(wallPos.reshape(1,len(wallPos),2), (len(position),1,1))
	wv = np.tile(wallVec.reshape(1,len(wallVec),2), (len(position),1,1))

	e = np.tile(position, (1,len(wallPos),1))
	ev = np.tile(velocity*timestep, (1,len(wallPos),1))
	
	t = np.cross(w - e, wv) / np.cross(ev, wv)
	u = np.cross(w - e, ev) / np.cross(ev, wv)

	mask = (t > 0) * (t < 1) * (u > 0) * (u < 1)

	wvn = np.tile(wallNorm.reshape(1,len(wallNorm),2), (len(velocity),1,1))

	if(np.sum(mask) >= 1) :
	
		reflections = ev / timestep - 2 * wvn * np.inner(velocity, wallNorm).reshape((len(velocity),len(wallNorm),1))
		maskedRefl = mask.reshape((len(velocity),len(wallNorm),1)) * reflections
		compressedRefl = np.sum(maskedRefl, axis=1,).reshape(len(velocity),1,2)

		invMask = np.sum(mask, axis=1,).reshape(len(velocity),1) < 1
		velocity[:,0] = invMask*velocity[:,0] + compressedRefl[:,0]


	position += velocity*timestep

	# Right periodic boundry
	rightMask = position[:,:,0] > width
	position[:,:,0] = (1-rightMask) * position[:,:,0] + rightMask * (position[:,:,0] - width)
	# Left periodic boundry
	leftMask = position[:,:,0] < 0
	position[:,:,0] = (1-leftMask) * position[:,:,0] + leftMask * (width + position[:,:,0])

	for i in range(0,min(len(position),MAX_PLOT_PATHS)) :
		if ~leftMask[i] and ~rightMask[i]:
			path[i][0].append(position[i][0][0])
			path[i][1].append(position[i][0][1])
		else :
			path[i][0].append(np.nan)
			path[i][1].append(np.nan)

def scatter(vtherm, velocities, timestep, meanTime) :
	scatterProb = 1 - np.exp(-1*timestep/meanTime)
	for velocity in velocities :
		if scatterProb > np.random.uniform(0,1) :
			angle = np.random.uniform(0,2*np.pi)
			velocity[0] = [vtherm*np.cos(angle),
				vtherm*np.sin(angle)]

def draw(plt, width, height, paths, edges) :
	ax = plt.gca()
	ax.set_xlim([0-width*0.5,width*1.5])
	ax.set_ylim([0-height*0.5,height*1.5])

	for path in paths :
		plt.plot(path[0], path[1])
	for edge in edges :
		#plt.arrow(edge[0], edge[1], edge[2], edge[3], head_width=2e-9, head_length=5e-9, width = 0.01e-9, fc='k', ec='k')
		plt.arrow(edge[0], edge[1], edge[2], edge[3], head_width=0.0, head_length=0.0, width = 0.01e-9, fc='k', ec='k')
	plt.show()


def heatMap(plt, width, height, positions) :
	lastX = []
	lastY = []
	for position in positions :
		# kinda a hack
		if position[0][1] < height and position[0][1] > 0:
			lastX.append(position[0][0])
			lastY.append(position[0][1])

	ax = plt.gca()
	ax.set_xlim([0,width])
	ax.set_ylim([0,height])

	ax.hist2d(lastX,lastY,bins=20)
	plt.show()

def trivialExample(numElectrons, numIterations):	
	start = time.clock()
	timestep = 5e-15

	width = 200e-9
	height = 100e-9
	vtherm = 1.87e5

	edges = [[0.0,height,width,0.0],[width,0.0,-width,0.0]]

	wallPosTemp = []
	wallVecTemp = []
	for edge in edges :
		wallPosTemp.append([edge[0:2]])
		wallVecTemp.append([edge[2:4]])
		wallPos = np.array(wallPosTemp)
		wallVec = np.array(wallVecTemp)

	wallNormTemp = []
	for wall in wallVec :
		length = np.sqrt(wall[0][0]**2 + wall[0][1]**2)
		x = wall[0][1] / length
		y = -1.0 * wall[0][0] / length
		wallNormTemp.append([[x,y]])

	wallNorm = np.array(wallNormTemp)	

	polys = []


	position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
	for i in range(numIterations) :
		iterate(plt, width, timestep, position, velocity, wallPos, wallVec, wallNorm, path)

	print(time.clock() - start)
	#draw(plt, width, height, path, edges)
	heatMap(plt, width, height, position)

def trapazoidExample(numElectrons, numIterations):	
	timestep = 5e-15

	width = 200e-9
	height = 100e-9
	vtherm = 1.87e5


	wallPos, wallVec, wallNorm, polys, edges = shapes.trapazoidBottleNeck(width,height)

	position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
	for i in range(numIterations) :
		scatter(vtherm, velocity, timestep, 0.2e-12)
		iterate(plt, width, timestep,  position, velocity, wallPos, wallVec, wallNorm, path)
	draw(plt, width, height, path, edges)
	heatMap(plt, width, height, position)	

def squareExample(numElectrons, numIterations):	
	timestep = 2e-15

	width = 200e-9
	height = 100e-9
	vtherm = 1.87e5


	wallPos, wallVec, wallNorm, polys, edges = shapes.squareBottleNeck(width,height)

	position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
	for i in range(numIterations) :
		scatter(vtherm, velocity, timestep, 0.2e-12)
		iterate(plt, width, timestep,  position, velocity, wallPos, wallVec, wallNorm, path)
		#print(i)
	draw(plt, width, height, path, edges)
	heatMap(plt, width, height, position)	
if __name__ == "__main__" :
	#trivialExample(1000,1000)
	squareExample(1000,1000)
