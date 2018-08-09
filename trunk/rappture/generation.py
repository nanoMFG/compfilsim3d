import numpy as np
import math
from distance import distance, simple_distance

"""Returns volume of polymer given width(of square side) and length."""
def get_polymer_volume(width_of_poly, length_of_poly):
    return (width_of_poly ** 2) * length_of_poly

"""Returns volume of cylindrical nanowire with diameter and length."""
def get_wire_volume(length, diameter):
    return math.pi * ((diameter / 2) ** 2) * length
  
"""Returns a random point within the domain. """ 
def get_initial_point(width_of_poly, length_of_poly):
    x = np.random.uniform(-(width_of_poly / 2), (width_of_poly / 2))
    y = np.random.uniform(0, length_of_poly)
    #y = np.random.uniform(0, 0)
    z = np.random.uniform(0, width_of_poly)
    return (x, y, z)

"""Returns an end point for the corresponding start_point within the domain. """
def get_end_point(p, mean, sd, phi_lower, phi_upper, theta_lower, theta_upper, 
                                                width_of_poly, length_of_poly):
    length = np.random.normal(mean, sd)
    phi = np.random.uniform(phi_lower, phi_upper)
    theta = verti_angle = np.random.uniform(theta_lower, theta_upper)
    x = p[0] + length * math.sin(theta) * math.cos(phi)
    y = p[1] + length * math.sin(theta) * math.sin(phi)
    z = p[2] + length * math.cos(theta)
    p_end = fix_domain(p, [x, y, z], -width_of_poly / 2, width_of_poly / 2, 0, length_of_poly, 0, width_of_poly)
    final_length = simple_distance(p, p_end)
    return tuple(p_end), final_length
  
"""Returns True if the line (point1, point2) has a distamce that is 
   greater than or equal to diameter away from all the lines defined by
   start_points, end_points."""  
def check_intersection(lines, point1, point2, diameter, int_tol):
    for line in lines:
        dist = distance(line, (point1, point2))
        if dist < diameter - int_tol:
            return False
    return True

"""Return True if the points (x, y, z) is not within the polymer. """
def check_domain(p):
    return p[1] > length_of_poly or p[2] > width_of_poly or p[2] < 0 \
                or abs(p[0]) > (width_of_poly / 2)

"""Given initial point p0 and a line (p0, p1), adjust p1 so it lies 
   within the domain x_lo, x_up, y_lo, y_up, z_lo, z_up."""
def fix_domain(p0, p1, x_lo, x_up, y_lo, y_up, z_lo, z_up):
    p0 = list(p0)
    p1 = list(p1)
    a = p1[0]-p0[0]
    b = p1[1]-p0[1]
    c = p1[2]-p0[2]
    if p1[0] > x_up:
        p1[0] = x_up
        t = (p1[0]-p0[0])/a
        p1[1] = p0[1] + t*b
        p1[2] = p0[2] + t*c
    elif p1[0] < x_lo:
        p1[0] = x_lo
        t = (p1[0]-p0[0])/a
        p1[1] = p0[1] + t*b
        p1[2] = p0[2] + t*c
    if p1[1] > y_up:
        p1[1] = y_up
        t = (p1[1]-p0[1])/b
        p1[0] = p0[0] + t*a
        p1[2] = p0[2] + t*c
    elif p1[1] < y_lo:
        p1[1] = y_lo
        t = (p1[1]-p0[1])/b
        p1[0] = p0[0] + t*a
        p1[2] = p0[2] + t*c  
    if p1[2] > z_up:
        p1[2] = z_up
        t = (p1[2]-p0[2])/c
        p1[0] = p0[0] + t*a
        p1[1] = p0[1] + t*b
    elif p1[2] < z_lo:
        p1[2] = z_lo
        t = (p1[2]-p0[2])/c
        p1[0] = p0[0] + t*a
        p1[1] = p0[1] + t*b
    return tuple(p1)
  
def get_cross_sectional_area(diameter):
    return (diameter / 2) ** 2 * math.pi
    
"""Returns a length=NUM_SUBSECTIONS list of list of lines where each line is [(x1, y1, z1), (x2, y2, z2)]."""
def generate(volume_fraction, width_of_poly, length_of_poly, diameter, 
            mean_nanowire_length, SD, theta_lower, 
            theta_upper, phi_lower, phi_upper, int_tol):
    volume_of_poly = get_polymer_volume(width_of_poly, length_of_poly)
    lines = []
    volume = 0
    while (volume / volume_of_poly) < volume_fraction:
        start = get_initial_point(width_of_poly, length_of_poly)
        end, length = get_end_point(start, mean_nanowire_length, SD, phi_lower, phi_upper, theta_lower, 
            theta_upper, width_of_poly, length_of_poly)
        while (int_tol != diameter and not check_intersection(lines, start, end, diameter, int_tol)):
            start = get_initial_point(width_of_poly, length_of_poly)
            end, length = get_end_point(start, mean_nanowire_length, SD, phi_lower, phi_upper, theta_lower, 
                theta_upper, width_of_poly, length_of_poly)
        # Enforce ordering of y values.
        if start[1] <= end[1]:
            lines.append([start, end])
        else:
            lines.append([end, start])
        volume += math.pi * ((diameter / 2) ** 2) * length
    return lines
    
def chop(lines, num_subsections, length_of_poly, width_of_poly):
    # Divide the lines along number of subsections.
    subsections = dict()
    for i in range(num_subsections):
        subsections[i] = []
    for lin in lines:
        y1 = lin[0][1]
        y2 = lin[1][1]
        assert(y1 <= y2)
        sublength = length_of_poly/num_subsections
        for i in range(num_subsections):
            cond = False
            if i == num_subsections - 1:
                cond = y2 >= sublength*i and y2 <= sublength*(i+1)
            else:
                cond = (y1 >= sublength*i and y1 < sublength*(i+1)) or \
               (y2 >= sublength*i and y2 < sublength*(i+1)) or \
               (y1 < sublength*i and y2 > sublength*(i+1))
            if cond:  
                p0 = lin[0]
                p1 = lin[1]
                p1 = fix_domain(p0, p1, -width_of_poly / 2, width_of_poly / 2, sublength*i, sublength*(i+1), 0, width_of_poly)
                p0 = fix_domain(p1, p0, -width_of_poly / 2, width_of_poly / 2, sublength*i, sublength*(i+1), 0, width_of_poly)
                subsections[i].append([p0, p1])
                #print("y values of %f and %f got placed in the %ith subsection" % (y1, y2, i))
    ret_val = []
    for i in subsections:
        ret_val.append(subsections[i])
    return ret_val