import matplotlib.pyplot as plt
import numpy as np
import shapes


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


	for i in range(numElectrons) :
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

def scatter(velocities, timestep, meanTime) :
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

if __name__ == "__main__" :
	numElectrons = 4
	timestep = 1

	width = 200e-9
	height = 100e-9
	vtherm = 0.5e-9

	wallPos, wallVec, wallNorm, polys, edges = shapes.trapazoidBottleNeck(width,height)

	position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
	for i in range(1000) :
		scatter(velocity, 0.1, 10)
		iterate(plt, position, velocity, wallPos, wallVec, wallNorm, path)
	draw(plt, width, height, path, edges)

