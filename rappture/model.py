import Rappture
import sys
import math
from generation import *
from graph import init_graph, find_connecting_cluster
from conductivity import total_conductivity
from spice import write_netlist

sys.stderr = open('model.err', 'w')
sys.stdout = open('model.out', 'w')
print("<0>")
io = Rappture.library(sys.argv[1])
# uncomment these for debugging

diameter = float(io.get('input.number(diameter).current'))*10**(-3)
mean_nanowire_length = float(io.get('input.number(length).current'))
volume_fraction = float(io.get('input.number(volume_fraction).current'))/100
# print("diameter is, ", diameter)
alignment_factor = float(io.get('input.number(alignment_factor).current'))

#diameter = (diameter)e-3
width_of_poly = 25
length_of_poly = 25
SD = .5
# theta_lower = 0
# theta_upper = 2 * math.pi
# phi_lower = 0
# phi_upper = 2 * math.pi
theta_lower = alignment_factor * -math.pi
theta_upper = alignment_factor * math.pi
phi_lower = alignment_factor * -math.pi
phi_upper = alignment_factor * math.pi
cross_sectional_area = (diameter / 2) ** 2 * math.pi
intersection_tolerance = diameter
num_subsections = 1




print("<1>")
# Uncomment this chunk to get percolation true/false outputted.
# lines = generate(volume_fraction, width_of_poly, 
#                                    length_of_poly, diameter, 
#                                    mean_nanowire_length, 
#                                    SD, theta_lower, theta_upper, 
#                                    phi_lower, phi_upper, 
#                                    intersection_tolerance)
# subsections = chop(lines, num_subsections, length_of_poly, width_of_poly) 
# H, lineswi = init_graph(lines, diameter*1) 
# segments = []
# for line in lines:
#     x = [line[0][0], line[1][0]]
#     y = [line[0][1], line[1][1]]
#     z = [line[0][2], line[1][2]]
#     segments.append([x, y, z])
# K, K_sources, K_sinks = find_connecting_cluster(segments, H, 0, length_of_poly)
# path = str(bool(K))
# io.put('output.string(path).current', path)

num_subsections_list = [3]
lineswi_dict, segments_dict, K_dict, K_sources_dict, K_sinks_dict = \
    repeat_generation(volume_fraction, width_of_poly, 
                                   length_of_poly, diameter, 
                                   mean_nanowire_length, 
                                   SD, theta_lower, theta_upper, 
                                   phi_lower, phi_upper, 
                                   intersection_tolerance, 10, num_subsections_list)
print("<2>")
lineswlen_dict = dict()
for num_subsections in num_subsections_list:
    i = num_subsections
    lineswlen_list = [process_intersection_points(lineswi) for lineswi in lineswi_dict[i]]
    lineswlen_dict[i] = lineswlen_list
# Using experiment (a) with Ag and R=0.3, p=0.4, C=1
print("<3>")
cond = total_conductivity(diameter*1e-6/2, diameter*1e-6, 0.3, diameter*1e-6, 1, 6.3e7, .4)


current_dict = dict()
for num_subsections in num_subsections_list:
    i = num_subsections
    print("Calling write_netlist for %i subsections" % (i))
    current_list = [write_netlist(K, lineswlen, K_sources, K_sinks, diameter*1e-6, segments, 5, j)
                   for K, lineswlen, K_sources, K_sinks, segments, j 
                   in zip(K_dict[i], lineswlen_dict[i], K_sources_dict[i], K_sinks_dict[i], segments_dict[i], range(i))] 
    current_dict[i] = current_list
print("<4>")
Rappture.result(io)
