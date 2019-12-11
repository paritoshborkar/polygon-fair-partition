import math
import numpy as np
from matplotlib import pyplot as plt

def main():

    numberOfCoords = 0
    origcoordinateList = []
    fileName = "datapoints.txt"
    with open(fileName) as file:
        for data in file:
            if data.__contains__(" "):
                data = data.split(" ")
                origcoordinateList.append((int(data[0].strip()), int(data[1].strip())))

    coordForTriang = sorted(origcoordinateList[:len(origcoordinateList)-1], key=lambda x: x[0])
    triangulatedGroups = []
    for index in range(len(coordForTriang)-2):
        triangulatedGroups.append(coordForTriang[index:index+3])

    L = []
    P1 = []
    P2 = origcoordinateList
    Per1 = 0
    Per2 = calculatePerimeter(origcoordinateList)

    for territory in range(len(triangulatedGroups)):
        region1Area = calculateArea(P1) + calculateArea(triangulatedGroups[territory])
        region2Area = calculateArea(P2) - calculateArea(triangulatedGroups[territory])

        if region1Area > region2Area:
            p = triangulatedGroups[territory][1]
            r = triangulatedGroups[territory][0]
            s = triangulatedGroups[territory][2]
            D1 = region1Area - region2Area
            tempQ = (D1 / (2 * calculateArea(triangulatedGroups[territory])))
            tempQ = [tempQ * (r[0] - s[0]), tempQ * (r[1] - s[1])]
            q = (tempQ[0] + s[0], tempQ[1] + s[1])
            L = [p, q]
            return [L, P1, P2, Per1, Per2]
        P1 = addTerritory(P1, triangulatedGroups[territory])
        P2 = removeTerritory(P2, triangulatedGroups[territory])
        if territory == 0:
            Per1 = Per1 + distance(triangulatedGroups[territory][1], triangulatedGroups[territory][2]) \
                   + distance(triangulatedGroups[territory][1], triangulatedGroups[territory][0])
            Per2 = Per2 - distance(triangulatedGroups[territory][1], triangulatedGroups[territory][2]) \
                   - distance(triangulatedGroups[territory][1], triangulatedGroups[territory][0])
        else:
            Per1 = Per1 + distance(triangulatedGroups[territory][1], triangulatedGroups[territory][2])
            Per2 = Per2 - distance(triangulatedGroups[territory][1], triangulatedGroups[territory][2])
    return None


def addTerritory(P=[], T=[]):
    if len(P) == 0:
        P = T
    else:
        P.pop()
        P.insert(P.index(T[1]), T[2])
    P.append(P[0])
    return P

def removeTerritory(P=[], T=[]):
    if len(P) > 3:
        P.remove(T[0])
        P[len(P)-1] = P[0]
        return P
    else:
        return []

def calculateArea(P=[]):
    if len(P) < 2:
        return 0
    area = 0
    for i in range(len(P)-1):
        point1 = P[i]
        point2 = P[i + 1]
        areaCal = (point2[0]-point1[0])*(point2[1]+point1[1])/2
        area -= areaCal
    point1 = P[-1]
    point2 = P[0]
    areaCal = (point2[0] - point1[0])*(point2[1] + point1[1])/2
    area -= areaCal
    return abs(area)

def calculatePerimeter(coorList=[]):
    perimeter = 0
    for i in range(len(coorList) - 1):
        point1 = coorList[i]
        point2 = coorList[i + 1]
        perimeter += distance(point1, point2)
    perimeter += distance(coorList[len(coorList) - 1], coorList[0])
    return perimeter

def distance(point1=[], point2=[]):
    return math.sqrt(math.pow(point2[1]-point1[1], 2) + math.pow(point2[0]-point1[0], 2))


result = main()
print(result)
x = []
y = []
resultP1 = []
resultP2 = []
P1 = result[1]
P1.insert(P1.index(result[0][0])+1, result[0][1])
P2 = result[2]
P2[0] = result[0][1]
P2[len(P2)-1] = result[0][1]

for p in P1:
    x.append(p[0])
    y.append(p[1])
plt.fill(x, y, color="b")
plt.plot(x, y)

x = []
y = []
for p in P2:
    x.append(p[0])
    y.append(p[1])
plt.fill(x, y, color="b")
plt.plot(x, y)

plt.show()
P1.reverse()
P2.reverse()
print(P1, P2)
print(calculateArea(P1), calculateArea(P2))