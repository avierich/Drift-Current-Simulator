import numpy as np

def trapazoidBottleNeck(width, height) :

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

	polys = [[[width/3,height,0.0,-height/4],
		[width/3,3*height/4,width/9,-height/8],
		[4*width/9,5*height/8,width/9,0.0],
		[5*width/9,5*height/8,width/9,height/8],
		[6*width/9,6*height/8,0.0,height/4]],
		[[2*width/3,0.0,0.0,height/4],
		[2*width/3,height/4,-width/9,height/8],
		[5*width/9,3*height/8,-width/9,0.0],
		[4*width/9,3*height/8,-width/9,-height/8],
		[3*width/9,2*height/8,0.0,-height/4]]]

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
	return wallPos, wallVec, wallNorm, polys, edges

def squareBottleNeck(width, height) :

	edges = [[0.0,height,2*width/5,0.0],
		[2*width/5,height,0.0,-2*height/5],
		[2*width/5,3*height/5,width/5,0],
		[3*width/5,3*height/5,0,2*height/5],
		[3*width/5,height,2*width/5,0],
		[width,0.0,-2*width/5,0],
		[3*width/5,0.0,0.0,2*height/5],
		[3*width/5,2*height/5,-width/5,0.0],
		[2*width/5,2*height/5,0.0,-2*height/5],
		[2*width/5,0.0,-2*width/5,0.0]]
			

	polys = [[[2*width/5,height,0.0,-2*height/5],
		[2*width/5,3*height/5,width/5,0],
		[3*width/5,3*height/5,0,2*height/5]],
		[[3*width/5,0.0,0.0,2*height/2],
		[3*width/5,2*height/5,-width/5,0.0],
		[2*width/5,2*height/5,0.0,-2*height/5]]]

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
	return wallPos, wallVec, wallNorm, polys, edges
