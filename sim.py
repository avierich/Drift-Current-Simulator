import matplotlib.pyplot as plt
import numpy as np


def initElectrons(numElectrons, width, height, vtherm) :
	tempPosition = []
	tempVelocity = []
	paths = []
	for i in range(numElectrons) :
		# Generate random position
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
	
	# Reflections
	#print(velocity.shape)
	#print(np.reshape(wallNorm,(1,len(wallNorm),2,1)))

	#print(np.dot(velocity,np.reshape(wallNorm,(1,len(wallNorm),2))))


	for i in range(0,len(t)) :
		for j in range(len(t[0])) :
			if(t[i][j] > 0 and t[i][j] < 1 and u[i][j] > 0 and u[i][j] < 1) :
				x = position[i][0][0] + velocity[i][0][0] * t[i][j]
				y = position[i][0][1] + velocity[i][0][1] * t[i][j]
				plt.plot([x,x], [y,y], 'x',markersize = 12)

	position += velocity
	for i in range(0,len(position)) :
		path[i][0].append(position[i][0][0])
		path[i][1].append(position[i][0][1])

def draw(plt, width, height, paths, edges) :
	ax = plt.gca()
	ax.set_xlim([0,width])
	ax.set_ylim([0,height])

	for path in paths :
		plt.plot(path[0], path[1])
	for edge in edges :
		plt.arrow(edge[0], edge[1], edge[2], edge[3], head_width=2e-9, head_length=5e-9, width = 0.01e-9, fc='k', ec='k')
	plt.show()

if __name__ == "__main__" :
	numElectrons = 30

	width = 200e-9
	height = 100e-9
	vtherm = 3e-9

	edges = [[150e-9,80e-9,0.0,-60e-9],[50e-9,80e-9,0.0,-60e-9],[125e-9,20e-9,-50e-9,60e-9]]

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

	position, velocity, path = initElectrons(numElectrons, width, height, vtherm)
	for i in range(1000) :
		iterate(plt, position, velocity, wallPos, wallVec, wallNorm, path)
	draw(plt, width, height, path, edges)
