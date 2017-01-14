import matplotlib.pyplot as plt
import numpy as np


def initElectrons(numElectrons, width, height, vtherm, polygons) :
	tempPosition = []
	tempVelocity = []
	paths = []
	for i in range(numElectrons) :
		# Generate random position
#		attemptPos = np.array([np.random.uniform(0, width),
#				np.random.uniform(0, height)]])
#		isValid = False
#		while ~isValid :
#			isValid = True
#			attemptPos = np.array([np.random.uniform(0, width),
#						np.random.uniform(0, height)]])
#			for polygon in polygons :
#				for wall in polygon :
#					originToPoint = np.array([attemptPos[0] - wall[0],attemptPos[1] - wall[1]])
#					wallVec = np.array([wall[2],wall[3]])
#					if np.dot(originToPoint, wallVec) > 0 :
#						isValid = False
					

		tempPosition.append([[np.random.uniform(0, width),
				np.random.uniform(0, height)]])

		angle = np.random.uniform(0,2*np.pi)
		tempVelocity.append([[vtherm*np.cos(angle),
				vtherm*np.sin(angle)]])


	for i in range(numElectrons) :
		paths.append([[],[]])

	return np.array(tempPosition), np.array(tempVelocity), paths

def iterate(plt,position, velocity, wallPos, wallVec, wallNorm, path) :

	# Collisions
	w = np.tile(wallPos.reshape(1,len(wallPos),2), (len(position),1,1))
	wv = np.tile(wallVec.reshape(1,len(wallVec),2), (len(position),1,1))

	e = np.tile(position, (1,len(wallPos),1))
	ev = np.tile(velocity, (1,len(wallPos),1))
	
	t = np.cross(w - e, wv) / np.cross(ev, wv)
	u = np.cross(w - e, ev) / np.cross(ev, wv)

	mask = (t > 0) * (t < 1) * (u > 0) * (u < 1)
	compressedMask = np.dot(mask, np.ones((len(mask[0]),1)))

	for i in range(0, len(mask)) :
		for j in range(0, len(mask[0])) :
			if(mask[i][j] == 1):
				velocity[i][0] = velocity[i][0] -  2*np.dot(velocity[i][0].reshape(1,2), wallNorm[j][0].reshape((2,1))) * wallNorm[j][0]

	for i in range(0,len(t)) :
		for j in range(len(t[0])) :
			if(t[i][j] > 0 and t[i][j] < 1 and u[i][j] > 0 and u[i][j] < 1) :
				x = position[i][0][0] + velocity[i][0][0] * t[i][j]
				y = position[i][0][1] + velocity[i][0][1] * t[i][j]
				#plt.plot([x,x], [y,y], 'x',markersize = 12)
				plt.plot([x,x], [y,y])

	position += velocity

	# Right periodic boundry
	rightMask = position[:,:,0] > width
	position[:,:,0] = (1-rightMask) * position[:,:,0] + rightMask * (position[:,:,0] - width)
	# Left periodic boundry
	leftMask = position[:,:,0] < 0
	position[:,:,0] = (1-leftMask) * position[:,:,0] + leftMask * (width + position[:,:,0])

	for i in range(0,len(position)) :
		if ~leftMask[i] and ~rightMask[i]:
			path[i][0].append(position[i][0][0])
			path[i][1].append(position[i][0][1])
		else :
			path[i][0].append(np.nan)
			path[i][1].append(np.nan)

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

if __name__ == "__main__" :
	numElectrons = 10

	width = 200e-9
	height = 100e-9
	vtherm = 0.5e-9

	edges = [[0.0,height,width/3,0.0],
		[width/3,height,0.0,-height/4],
		[width/3,3*height/4,width/9,-height/8],
		[4*width/9,5*height/8,width/9,0.0],
		[5*width/9,5*height/8,width/9,height/8],
		[6*width/9,6*height/8,0.0,height/4],
		[6*width/9,height,width/3,0.0],
		[width,0.0,-width/3,0.0],
		[2*width/3,0.0,0.0,height/4],
		[2*width/3,height/4,-width/9,height/8],
		[5*width/9,3*height/8,-width/9,0.0],
		[4*width/9,3*height/8,-width/9,-height/8],
		[3*width/9,2*height/8,0.0,-height/4],
		[3*width/9,0.0,-width/3,0.0]]

	polys = [[[0.0,height,width/3,0.0],
		[width/3,height,0.0,-height/4],
		[width/3,3*height/4,width/9,-height/8],
		[4*width/9,5*height/8,width/9,0.0],
		[5*width/9,5*height/8,width/9,height/8],
		[6*width/9,6*height/8,0.0,height/4],
		[6*width/9,height,width/3,0.0]],
		[[width,0.0,-width/3,0.0],
		[2*width/3,0.0,0.0,height/4],
		[2*width/3,height/4,-width/9,height/8],
		[5*width/9,3*height/8,-width/9,0.0],
		[4*width/9,3*height/8,-width/9,-height/8],
		[3*width/9,2*height/8,0.0,-height/4],
		[3*width/9,0.0,-width/3,0.0]]]

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

	position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
	for i in range(1000) :
		iterate(plt, position, velocity, wallPos, wallVec, wallNorm, path)
	draw(plt, width, height, path, edges)
