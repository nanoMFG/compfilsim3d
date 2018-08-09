import Rappture
import sys
import math
from generation import generate, chop
from graph import init_graph

io = Rappture.library(sys.argv[1])

diameter = float(io.get('input.number(diameter).current'))
mean_nanowire_length = float(io.get('input.number(length).current'))
volume_fraction = float(io.get('input.number(volume_fraction).current'))

width_of_poly = 25
length_of_poly = 25
SD = .5
theta_lower = 0
theta_upper = 2 * math.pi
phi_lower = 0
phi_upper = 2 * math.pi
cross_sectional_area = (diameter / 2) ** 2 * math.pi
intersection_tolerance = diameter
num_subsections = 1
lines = generate(volume_fraction, width_of_poly, 
                                   length_of_poly, diameter, 
                                   mean_nanowire_length, 
                                   SD, theta_lower, theta_upper, 
                                   phi_lower, phi_upper, 
                                   intersection_tolerance)
subsections = chop(lines, num_subsections, length_of_poly, width_of_poly) 

H, lineswi = init_graph(lines, diameter*1) 
for line in lines:
    x = [line[0][0], line[1][0]]
    y = [line[0][1], line[1][1]]
    z = [line[0][2], line[1][2]]
    segments.append([x, y, z])
K, K_sources, K_sinks = find_connecting_cluster(segments, H, 0, length_of_poly)

path = bool(K)
io.put('output.string(path).current')
# print 'formula = %s' % formula
# npts = 100

# io.put('output.curve(result).about.label','Formula: Y vs X',append=0)
# io.put('output.curve(result).yaxis.label','Y')
# io.put('output.curve(result).xaxis.label','X')

# for i in range(npts):
#     x = (xmax-xmin)/npts * i + xmin;
#     y = eval(formula)
#     io.put('output.curve(result).component.xy', '%g %g\n' % (x,y), append=1)

Rappture.result(io)