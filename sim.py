import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import shapes
import time

T = 300
m = 0.26 * 9.11e-31
k = 1.38e-23

MAX_PLOT_PATHS = 20

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
        

        tempVelocity.append([[np.random.randn()*np.sqrt(k*T/m),
                                np.random.randn()*np.sqrt(k*T/m)]])


    for i in range(min(numElectrons,MAX_PLOT_PATHS)) :
        paths.append([[],[]])

    return np.array(tempPosition), np.array(tempVelocity), paths

# This function checks if a given point lies in any of the given polygons
def isPointInPolys(polygons,point):
    for polygon in polygons :
        numIntersections = 0
        for edge in polygon :
            wallPos = np.array([edge[0], edge[1]])
            wallVec = np.array([edge[2], edge[3]])
            pointPos = np.array([-100e-9,-100e-9])
            pointVec = np.array([point[0] + 100e-9, point[1] +100e-9])

            t = np.cross(wallPos - pointPos, wallVec) / np.cross(pointVec, wallVec)
            u = np.cross(wallPos - pointPos, pointVec) / np.cross(pointVec, wallVec)
            
            numIntersections += t > 0 and t < 1 and u > 0 and u < 1

        if numIntersections % 2 :
            return True
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
            velocity[0] = [np.random.randn()*np.sqrt(k*T/m),
                np.random.randn()*np.sqrt(k*T/m)]

def draw(plt, width, height, paths, edges) :
    ax = plt.gca()
    ax.set_xlim([0,width])
    ax.set_ylim([0,height])

    for path in paths :
        plt.plot(path[0], path[1])
    for edge in edges :
        plt.arrow(edge[0], edge[1], edge[2], edge[3], head_width=0.0, head_length=0.0, width = 0.01e-9, fc='k', ec='k')
    plt.title("Electron Paths")
    plt.xlabel("X Position (m)")
    plt.ylabel("Y Position (m)")
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

    plt.hist2d(lastX,lastY,bins=20)
    plt.colorbar(norm=mcolors.NoNorm)
    plt.title("Final Electron Position")
    plt.xlabel("X Position (m)")
    plt.ylabel("Y Position (m)")
    plt.show()

def velocityHistogram(plt, velocities) :
    totalVelocity = []
    for velocity in velocities :
        totalVelocity.append((velocity[0][0]**2+velocity[0][1]**2)**0.5)

    plt.hist(totalVelocity, 10)
    plt.title("Final Electron Velocity")
    plt.xlabel("Electron Velocity (m/s)")
    plt.ylabel("Number of Electrons")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.show()

def plotTemperature(plt, temperature):
    plt.plot(temperature)
    plt.title("Semiconductor Temperature")
    plt.xlabel("Iteration Number")
    plt.ylabel("Temperature (K)")
    plt.show()

def measureTemperature(temperatures, velocities):
    temperature = 0
    for velocity in velocities :
        temperature += velocity[0][0]**2+velocity[0][1]**2

    temperatures.append(temperature*m/(k*2*len(velocities)))

def simulate(numElectrons, numIterations, world, scatterElectrons = True, showProgress = False) :
    timestep = 2e-15
    width = 200e-9
    height = 100e-9
    vtherm = np.sqrt((2*k*T)/m)

    temperature = []

    wallPos, wallVec, wallNorm, polys, edges = world(width,height)

    position, velocity, path = initElectrons(numElectrons, width, height, vtherm, polys)
    for i in range(numIterations) :
        if scatterElectrons :
            scatter(vtherm, velocity, timestep, 1.2e-12)
        iterate(plt, width, timestep,  position, velocity, wallPos, wallVec, wallNorm, path)
        measureTemperature(temperature, velocity)
        if showProgress :
            print(i)

    draw(plt, width, height, path, edges)
    heatMap(plt, width, height, position)
    velocityHistogram(plt, velocity)
    plotTemperature(plt, temperature)

if __name__ == "__main__" :
    simulate(1000,1000,shapes.squareBottleNeck,showProgress = True)

