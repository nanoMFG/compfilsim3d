from conductivity import total_conductivity
import numpy as np

def calc_resistance(total_cond, cross_sectional_area, length):
    resistivity = 1 / total_cond
    resistance = (resistivity * length) / cross_sectional_area
    return resistance

def calc_resistvity(R, cross_sectional_area, length):
    return R*cross_sectional_area/length

# Compute the length of the segment [(x1, x2), (y1, y2), (z1, z2)].
def get_length(seg):
    x1 = seg[0][0]
    x2 = seg[0][1]
    y1 = seg[1][0]
    y2 = seg[1][1]
    z1 = seg[2][0]
    z2 = seg[2][1] 
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    

# Given the connecting cluster K, runs pyspice on the graph. LINESWLEN is a dict 
# {line#: [d1, (p1, n1), d2, (p2, n2), d3]} where p# are intersecting points and 
# n# is the line number.
def write_netlist(K, lineswlen, source_nodes, sink_nodes, diameter, segments, source_voltage, subsection_num):
    #print("sink nodes:", sink_nodes)
    # First, create a dict where each element represents a wire and contains its
    # nodes and subdivided lengths. For example, value = [1a, d1, {1, 123}, d2, sink].
    # Endpoints are either #a, #b, source, or gnd. Set pairs are unordered.
    # Inner nodes are size 2 of sets, representing which wires are intersecting.
    wireswnodes = dict()
    for node in K:
        wirevalue = []
        if node in source_nodes:
            wirevalue.append('source')
        else:
            wirevalue.append(str(node) + 'a')
        wiredata = lineswlen[node]
        for item in wiredata:
            if not isinstance(item, tuple):
                wirevalue.append(item)
            else:
                if node < item[1]:
                    wirevalue.append(str(node) + "-" + str(item[1]))
                else:
                    wirevalue.append(str(item[1])+ "-" + str(node))
        if node in sink_nodes:
            wirevalue.append('gnd')
        else:
            wirevalue.append(str(node) + 'b')
        assert(len(wirevalue) >= 3 and len(wirevalue)%2 == 1) #Sanity check
        wireswnodes[node] = wirevalue
        
    # Next, create the circuit.
    Nodes = {}
    Wires = {}
    source_name = "Test" + str(subsection_num) + ".cir"
    fh = open(source_name, "w")
    fh.write("***Nanowire Composite Simulation w/ Resistors Only\n")
    fh.write("Vinput source 0 %dV\n" % (source_voltage))     
    cond = total_conductivity(diameter*1e-6/2, diameter*1e-6, 0.3, diameter*1e-6, 1, 6.3e7, .4)
    for wire, value in wireswnodes.items():
        for i in range(0, len(value)-1, 2):
            r_name = str(wire)+"_"+str(int(i/2)) # Name resistors w#_s#, w# is wire and s# is section number of the wire.
            n1 = value[i]
            sublength = value[i+1]*1e-6
            subr = calc_resistance(cond, np.pi*(diameter/2)**2, sublength)
            n2 = value[i+2]
            if n2 == 'gnd' and (subr != 0):
                Wires[r_name] = [r_name, n1, 0, subr]
            else:
                if subr != 0:
                    Wires[r_name] = [r_name, n1, n2, subr]
            if n1 in Nodes and subr != 0:
                update = Nodes[n1]
                update.append(r_name)
            else:
                if subr != 0:
                    Nodes[n1] = [r_name]
            if n2 in Nodes and subr != 0:
                update = Nodes[n2]
                update.append(r_name)
            else:
                if subr != 0:
                    Nodes[n2] = [r_name]
    deletewires = []
    for node, wire_list in Nodes.items():
        if len(wire_list) == 1:
            wire_name = wire_list[0]
            wire = Wires[wire_name]
            if wire[1] == node:
                compare = wire[2]
            else: 
                compare = wire[1]
            if compare != 0 and len(Nodes[compare]) == 1:
                if wire_name not in deletewires:
                    deletewires.append(wire_name)
    
    for wire in deletewires:
        print("deleting ", wire)
        del Wires[wire]
                
    for wire, wire_list in Wires.items():
        if wire_list[3] != 0:         
            fh.write("R%s %s %s %E\n" % (wire_list[0], wire_list[1], wire_list[2], wire_list[3]))

    fh.close()