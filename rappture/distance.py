import numpy as np
import math

"""Distance returns the 3D distance between line1=(start1, end1) and
line2=(start2, end2)."""
def distance(line1, line2):
    P = np.array(line1[0])
    Q = np.array(line1[1])
    R = np.array(line2[0])
    S = np.array(line2[1])
    P1mP0 = Q - P
    Q1mQ0 = S - R
    P0mQ0 = P - R
    a = np.dot(P1mP0, P1mP0)
    b = np.dot(P1mP0, Q1mQ0)
    c = np.dot(Q1mQ0, Q1mQ0)         
    d = np.dot(P1mP0, P0mQ0)
    e= np.dot(Q1mQ0, P0mQ0)
    det = a * c - b * b
    if det > 0:
        bte = b * e
        ctd = c * d
        if bte <= ctd: # s <= 0
            s = 0
            if e <= 0:  # t <= 0
                t = 0
                nd = -d
                if nd >= a:
                    s = 1
                elif nd > 0:
                    s = nd / a
            elif e < c:  #0 < t < 1 and region 5
                t = e / c
            else: # t >= 1
                t = 1   # region 4
                bmd = b - d
                if bmd >= a:
                    s = 1
                elif bmd > 0:
                    s = bmd / a
        else: 
            s = bte - ctd
            if s >= det:
                s = 1
                bpe = b + e
                if bpe <= 0: # region 8
                    t = 0
                    nd = -d
                    if nd <= 0:
                        s = 0
                    elif nd < a:
                        s = nd / a
                elif bpe < c: # region 1
                    t = bpe / c
                else: # region 2
                    t = 1
                    bmd = b- d
                    if bmd <= 0:
                        s = 0
                    elif bmd < a:
                        s = bmd / a
            else:
                ate = a * e
                btd = b * d
                if ate <= btd: # region 7 
                    t = 0
                    nd = -d
                    if nd <= 0:
                        s = 0
                    elif nd >= a:
                        s = 1
                    else:
                        s = nd / a
                else: 
                    t = ate - btd
                    if t >= det: # region 3
                        t = 1
                        bmd = b - d
                        if bmd <= 0:
                            s = 0
                        elif bmd >= a:
                            s = 1
                        else:
                            s = bmd / a
                    else:  # region 0
                        s /= det
                        t /= det
    else: 
        if e <= 0:
            t = 0
            nd = -d
            if nd <= 0: # region 6
                s = 0
            elif nd >=a: # region 8
                s = 1
            else: # region 7
                s = nd / a
        elif e >= c:
            t = 1
            bmd = b - d
            if bmd <= 0: # reigon 4
                s = 0
            elif bmd >= a: # region 2
                s = 1
            else:  # region 3
                s = bmd / a
        else:       
            """The point (0,e/c) is on the line and domain, so we have one
            point at which R is a minimum."""
            s = 0
            t = e / c

    result_parameter = [s, t]
    result_closest = [P + s * P1mP0, R + t * Q1mQ0]
    diff = result_closest[1] - result_closest[0]
    dist = np.dot(diff, diff)
    dist = math.sqrt(dist)
  
    return dist

# Compute the euclidian distance between p1=(x1,y1,z1) and p2=(x2,y2,z2)
def simple_distance(p1, p2):
    x1 = p1[0]
    x2 = p2[0]
    y1 = p1[1]
    y2 = p2[1]
    z1 = p1[2]
    z2 = p2[2]
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def clamppoint(point, line):
    x = np.array(line[0])
    y = np.array(line[1])
    point = np.array(point)
    clamped_point = [0, 0, 0]
    if x[0] <= y[0]:
        minX = x[0]
        maxX = y[0]
    else:
        minX = y[0]
        maxX = x[0]
    if x[1] <= y[1]:
        minY = x[1]
        maxY = y[1]
    else:
        minY = y[1]
        maxY = x[1]    
    if x[2] <= y[2]:
        minZ = x[2]
        maxZ = y[2]
    else:
        minZ = y[2]
        maxZ = x[2]
    if point[0] < minX:
        clamped_point[0] = minX
    else:
        if point[0] > maxX:
            clamped_point[0] = maxX
        else:
            clamped_point[0] = point[0]
    if point[1] < minY:
        clamped_point[1] = minY
    else:
        if point[1] > maxY:
            clamped_point[1] = maxY
        else:
            clamped_point[1] = point[1]
    if point[2] < minZ:
        clamped_point[2] = minZ
    else:
        if point[2] > maxZ:
            clamped_point[2] = maxZ
        else:
            clamped_point[2] = point[2]
    return tuple(clamped_point)

def distBetweenLines(line1, line2):
    p1 = np.array(line1[0])
    p2 = np.array(line1[1])
    p3 = np.array(line2[0])
    p4 = np.array(line2[1])
    d1 = p2 - p1
    d2 = p4 - p3
    eq1nCoeff = np.dot(d1, d2)
    eq1mCoeff = (-(d1[0] ** 2)) - (d1[1] ** 2) - (d1[2] ** 2)
    eq1Const = ((d1[0] * p3[0]) - (d1[0] * p1[0]) + (d1[1] * p3[1]) - (d1[1] * p1[1]) + (d1[2] * p3[2]) - (d1[2] * p1[2]))
    eq2nCoeff = (d2[0] ** 2) + (d2[1] ** 2) + (d2[2] ** 2)
    eq2mCoeff = -(d1[0] * d2[0]) - (d1[1] * d2[1]) - (d1[2] * d2[2]);
    eq2Const = ((d2[0] * p3[0]) - (d2[0] * p1[0]) + (d2[1] * p3[1]) - (d2[1] * p1[1]) + (d2[2] * p3[2]) - (d2[2] * p1[2]))
    M = [[ eq1nCoeff, eq1mCoeff, -eq1Const ], [ eq2nCoeff, eq2mCoeff, -eq2Const ]]
    a = np.array([[eq1nCoeff, eq1mCoeff], [eq2nCoeff, eq2mCoeff]])
    #for parallel lines
    if np.linalg.det(a) == 0:
        return [p1, p3]
    b = np.array([-eq1Const, -eq2Const])
    x = np.linalg.solve(a,b)
    n = x[0]
    m = x[1]
    i1 = [p1[0] + (m * d1[0]), p1[1] + (m * d1[1]), p1[2] + (m * d1[2])]
    i2 = [p3[0] + (n * d2[0]), p3[1] + (n * d2[1]), p3[2] + (n * d2[2])]
    i1clamped = clamppoint(i1, line1)
    i2clamped = clamppoint(i2, line2)
    return (i1clamped, i2clamped)
