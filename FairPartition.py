import sys
import math
import numpy as np
import math
from matplotlib import pyplot as plt


def getVertices(filename):
    vertices = []
    with open(filename) as file:
        for data in file:
            if data.__contains__(" "):
                data = data.split(" ")
                vertices.append([int(data[0].strip()), int(data[1].strip())])

    return vertices


def main(vertices):
    coordForTriang = sorted(vertices[:len(vertices) - 1], key=lambda x: x[0])
    triangulatedGroups = []
    for index in range(len(coordForTriang) - 2):
        triangulatedGroups.append(coordForTriang[index:index + 3])

    L = []
    P1 = []
    P2 = [x for x in vertices]
    Per1 = 0
    Per2 = calculatePerimeter(vertices)

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
            q = [tempQ[0] + s[0], tempQ[1] + s[1]]
            L = [p, q, s, r]
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
        P[len(P) - 1] = P[0]
        return P
    else:
        return []


def calculateArea(P=[]):
    if len(P) < 2:
        return 0
    area = 0
    for i in range(len(P) - 1):
        point1 = P[i]
        point2 = P[i + 1]
        areaCal = (point2[0] - point1[0]) * (point2[1] + point1[1]) / 2
        area -= areaCal
    point1 = P[-1]
    point2 = P[0]
    areaCal = (point2[0] - point1[0]) * (point2[1] + point1[1]) / 2
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
    return math.sqrt(math.pow(point2[1] - point1[1], 2) + math.pow(point2[0] - point1[0], 2))


"""
************************************************************************************
"""


class AB:
    __slots__ = 'p', 'q', 'r', 's'

    def __init__(self, p, q, r, s):
        self.p = p
        self.q = q
        self.r = r
        self.s = s


def area_polygon(vertices):
    result = 0
    for i in range(len(vertices)):
        current_vertex = vertices[i]
        next_vertex = vertices[(i + 1) % len(vertices)]

        result += ((next_vertex[0] - current_vertex[0]) * (current_vertex[1] + next_vertex[1])) / 2

    return abs(result)


def perimeter_polygon(vertices):
    result = 0
    for i in range(len(vertices)):
        current_vertex = vertices[i]
        next_vertex = vertices[(i + 1) % len(vertices)]

        result += segment_length(current_vertex, next_vertex)

    return result


def segment_length(vertex1, vertex2):
    return math.sqrt((vertex2[0] - vertex1[0]) ** 2 + (vertex2[1] - vertex1[1]) ** 2)


def cross_product(vector1, vector2):
    return (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])


def sin_vertices(vertice1, vertice2, vertice3):
    length1 = segment_length(vertice1, vertice2)
    length2 = segment_length(vertice3, vertice2)

    crossproduct = cross_product([vertice1[0] - vertice2[0], vertice1[1] - vertice2[1]],
                                 [vertice3[0] - vertice2[0], vertice3[1] - vertice2[1]])
    return crossproduct / (length1 * length2)


def intersection(segment1, segment2):
    a1 = (segment1[1][1] - segment1[0][1]) / (segment1[1][0] - segment1[0][0])
    a2 = (segment2[1][1] - segment2[0][1]) / (segment2[1][0] - segment2[0][0])

    b1 = (a1 * segment1[0][0]) - segment1[0][1]
    b2 = (a2 * segment2[0][0]) - segment2[0][1]

    x = (b2 - b1) / (a2 - a1)
    y = (a1 * x) - b1

    return [x, y]


def update(vertices, area_bisector, set_bisectors):
    p = area_bisector.p
    q = area_bisector.q
    r = area_bisector.r
    s = area_bisector.s

    p_index = vertices.index(p)
    p1 = vertices[(p_index - 1) % len(vertices)]
    p2 = vertices[(p_index + 1) % len(vertices)]

    o1 = intersection([p, q], [s, p1])  # pq intersection sp1
    o2 = intersection([p, q], [r, p2])  # pq intersection rp2

    p1_prime = None
    p2_prime = None
    q1_prime = None
    q2_prime = None

    L1 = None
    L2 = None
    L3 = None
    L4 = None
    if area_polygon([q, r, o2]) <= area_polygon([p, p2, o2]):  # pq
        num = ((area_polygon([p, p2, o2]) - area_polygon([q, r, o2])) /
               cross_product([r[0] - p2[0], r[1] - p2[1]], [p[0] - p2[0], p[1] - p2[1]]))

        vector = [p[0] - p2[0], p[1] - p2[1]]
        offset_vector = [num * element for element in vector]

        p1_prime = [p2[0] + offset_vector[0], p2[1] + offset_vector[1]]
        L1 = AB(r, p1_prime, p2, p)
        # L1 = [r, p1_prime]
    if area_polygon([q, s, o1]) <= area_polygon([p, p1, o1]):
        num = ((area_polygon([p, p1, o1]) - area_polygon([q, s, o1])) /
               cross_product([p[0] - p1[0], p[1] - p1[1]], [s[0] - p1[0], s[1] - p1[1]]))

        vector = [p[0] - p1[0], p[1] - p1[1]]
        offset_vector = [num * element for element in vector]

        p2_prime = [p1[0] + offset_vector[0], p1[1] + offset_vector[1]]
        L2 = AB(s, p2_prime, p, p1)
        # L2 = [s, p2_prime]
    if area_polygon([q, s, o1]) >= area_polygon([p, p1, o1]):
        sq = segment_length(s, q)
        sp1 = segment_length(s, p1)
        sin_rsp1 = sin_vertices(r, s, p1)
        num = (area_polygon([q, s, o1]) - (area_polygon([p, p1, o1])) / sq * sp1 * sin_rsp1)

        vector = [q[0] - s[0], q[1] - s[1]]
        offset_vector = [num * element for element in vector]

        q1_prime = [s[0] + offset_vector[0], s[1] + offset_vector[1]]
        L3 = AB(p1, q1_prime, r, s)
        # L3 = [p1, q1_prime]
    if area_polygon([q, r, o2]) >= area_polygon([p, p2, o2]):
        rq = segment_length(r, q)
        rp2 = segment_length(r, p2)
        sin_srp2 = sin_vertices(s, r, p2)

        num = (area_polygon([q, r, o2]) - (area_polygon([p, p2, o2])) / rq * rp2 * sin_srp2)

        vector = [q[0] - r[0], q[1] - r[1]]
        offset_vector = [num * element for element in vector]

        q2_prime = [r[0] + offset_vector[0], r[1] + offset_vector[1]]
        L4 = AB(p2, q2_prime, r, s)
        # L4 = [p2, q2_prime]

    L1_prime = None
    L2_prime = None

    if L1 is not None:
        L1_prime = L1
    else:
        L1_prime = L4

    if L2 is not None:
        L2_prime = L2
    else:
        L2_prime = L3

    # if p1_prime is None:
    #     p1_prime = q2_prime # q1_prime in paper
    #     r = p2
    # if p2_prime is None:
    #     p2_prime = q1_prime # q2_prime in paper
    #     s = p1

    a = None
    b = None

    P1 = [q]
    P2 = []
    for i in range(0, (vertices.index(p) - vertices.index(r) + 1 + len(vertices)) % len(vertices)):
        P1.append(vertices[(vertices.index(r) + i) % len(vertices)])
    for i in range(0, (vertices.index(s) - vertices.index(p) + 1 + len(vertices)) % len(vertices)):
        P2.append(vertices[(vertices.index(p) + i) % len(vertices)])
    P2.append(q)

    P11 = [L1_prime.q]
    P12 = []
    for i in range(0, (vertices.index(L1_prime.p) - vertices.index(L1_prime.r) + 1 + len(vertices)) % len(vertices)):
        P11.append(vertices[(vertices.index(L1_prime.r) + i) % len(vertices)])
    for i in range(0, (vertices.index(L1_prime.s) - vertices.index(L1_prime.p) + 1 + len(vertices)) % len(vertices)):
        P12.append(vertices[(vertices.index(L1_prime.p) + i) % len(vertices)])
    P12.append(L1_prime.q)

    P21 = [L2_prime.q]
    P22 = []
    for i in range(0, (vertices.index(L2_prime.p) - vertices.index(L2_prime.r) + 1 + len(vertices)) % len(vertices)):
        P21.append(vertices[(vertices.index(L2_prime.r) + i) % len(vertices)])
    for i in range(0, (vertices.index(L2_prime.s) - vertices.index(L2_prime.p) + 1 + len(vertices)) % len(vertices)):
        P22.append(vertices[(vertices.index(L2_prime.p) + i) % len(vertices)])
    P22.append(L2_prime.q)

    return [P11, P12, P21, P22]


"""
************************************************************************************
"""
def plotPts(P1=[], P2=[]):
    x = []
    y = []
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


filename = sys.argv[1]
vertices = getVertices(filename)
result = main(vertices)
L = result[0]
L = AB(L[0], L[1], L[2], L[3])
vertices.pop()

P1 = result[1]
P1.insert(P1.index(result[0][0])+1, result[0][1])
P2 = result[2]
P2[0] = result[0][1]
P2[len(P2)-1] = result[0][1]
plotPts(P1, P2)
resultPts = update(vertices, L, [])
resultPts[0].append(resultPts[0][0])
resultPts[1].append(resultPts[1][0])
resultPts[2].append(resultPts[2][0])
resultPts[3].append(resultPts[3][0])
print(resultPts)
plotPts(resultPts[0], resultPts[1])
plotPts(resultPts[2], resultPts[3])
# for _ in range(len(vertices)):
#     L, P1, P2, isDone = update(vertices, L, [])
#     print(isDone)
#     if isDone:
#         break
# print(L.p, ' ', L.q)
# print(P1)
# print(P2)
# print(perimeter_polygon(P1), perimeter_polygon(P2))
# print(calculateArea(P1), calculateArea(P2))

# Code to calculate exact points
# 	dp = perimeter_polygon(P1) - perimeter_polygon(P2)
#     dp1 = perimeter_polygon(P11) - perimeter_polygon(P12)
#     dp2 = perimeter_polygon(P21) - perimeter_polygon(P22)
#
#     L_prime = None
#     P1_prime = None
#     P2_prime = None
#
#     if dp == 0:
#         a = p
#         b = q
#         L_prime = AB(p, q, r, s)
#         P1_prime = P1
#         P2_prime = P2
#
#         return L_prime, P1_prime, P2_prime, True
#     elif dp * dp1 < 0 or dp * dp2 < 0:
#         d = None
#         c = None
#         if dp * dp1 < 0:  # dp = dp1 in paper
#             d = L1_prime.p
#             c = L1_prime.q
#         else:  # dp = dp2 in paper
#             d = L2_prime.p
#             c = L2_prime.q
#
#         area_dpq = area_polygon([d, p, q])
#         area_cpq = area_polygon([c, p, q])
#         A = area_dpq - area_cpq
#         B = (area_cpq * (dp/2 + segment_length(q, d))) - (area_dpq * (dp/2 + segment_length(p, c)))
#         C = - area_cpq * dp * segment_length(q, d)
#
#         qb1 = None
#         qb2 = None
#         qb = None
#         if area_dpq == area_cpq:
#             qb1 = -C / B
#             qb2 = qb1
#             qb = -C / B
#         else:
#             qb1 = (-B + math.sqrt(abs((B ** 2) - (4 * A * C)))) / (2 * A)
#             qb2 = (-B - math.sqrt(abs((B ** 2) - (4 * A * C)))) / (2 * A)
#             qb = (-B + math.sqrt(abs((B ** 2) - (4 * A * C)))) / (2 * A)
#
#         qb = qb1
#         b = [q[i] + (d[i] - q[i]) * qb for i in range(len(q))]
#         pa = qb - dp/2  # dp/2 in code part of paper
#         a = [p[i] + (c[i] - p[i]) * pa for i in range(len(p))]
#         O = intersection([a, b], [p, q])  # intersection of ab and pq
#
#         if abs(area_polygon([a, O, p]) - area_polygon([b, O, q])) > 0.1:
#             qb = qb2
#             b = [q[i] + (d[i] - q[i]) * qb for i in range(len(q))]
#             pa = qb - dp/2  # dp/2 in code part of paper
#             a = [p[i] + (c[i] - p[i]) * pa for i in range(len(p))]
#
#         L_prime = AB(a, b, q, d)
#         P1_prime = [b] + P1
#         P2_prime = [a] + P2
#         P1_prime.pop()
#         P2_prime.pop()
#
#         return L_prime, P1_prime, P2_prime, True
#     else:
#         if dp1 < dp2:
#             L_prime = L1_prime
#             P1_prime = P11
#             P2_prime = P12
#         else:
#             L_prime = L2_prime
#             P1_prime = P21
#             P2_prime = P22
#
#         return L_prime, P1_prime, P2_prime, False