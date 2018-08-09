import networkx as nx

""" Returns a networkx graph G formed from V = the lines in format 
 [(x1, y1, z1), (x2, y2, z2)] and E = distances between the lines
 that are within tolerance. Also returns a modified ordered lines with 
 format [endpoint1, endpoint2, (p1, n1), ...] where p1 is an 
 intersection point with line n.
"""
def init_graph(lines, tolerance):
    G = nx.Graph()
    G.add_nodes_from(range(len(lines)))
    lineswi = list(lines) # Lines with intersection
    for i in range(len(lines)):
        for j in range(i):
            l1 = lines[i]
            l2 = lines[j]
            d = distance(l1, l2)
            if d <= tolerance:
                G.add_edge(i, j, weight=d*1e-6)  # Units unscaled
                p = distBetweenLines(l1, l2)
                lines[i].append((p[0], j))
                lines[j].append((p[1], i))
    return G, lineswi

def is_between(b1, x, b2):
     return (b1 <= x and x <= b2) or (b2 <= x and x <= b1)

""" Given the length of the cross-section, a list of segments 
 [(x1, x2), (y1, y2), (z1, z2)], and graph G return a 
 connecting cluster from source(y~0) to sink(y~length). If none, 
 return None. 
"""
def find_connecting_cluster(segs, G, source_plane, sink_plane):
  
    n = range(len(segs))
    source_segs = [i for i in n if segs[i][1][0] == source_plane or 
                                   segs[i][1][1] == source_plane]
    sink_segs = [i for i in n if segs[i][1][0] == sink_plane or 
                                 segs[i][1][1] == sink_plane]
  
    if not source_segs or not sink_segs:
        print("No nanotubes intersecting source and/or sink planes.")
        return None, None, None
  
    # Iteratively grow the connecting cluster beginning from the source.
    H = nx.Graph()
    H.add_nodes_from(source_segs)
    new_segs = source_segs
    while new_segs:
        temp = list(new_segs)
        new_segs = []
        for seg in temp:
            for v in G.neighbors(seg):
                if v not in sink_segs and v not in H.nodes():
                    new_segs.append(v)
                H.add_edge(seg, v, weight=G[seg][v]['weight'])  
  
    # Iteratively grow from the sink to get the final cluster.
    K = nx.Graph()
    K_sinks, K_sources = set(), set()
    new_segs = [v for v in H.nodes() if v in sink_segs]
    K_sinks = set(new_segs)
    K.add_nodes_from(new_segs)
    path_exists = False
    while new_segs:
        temp = list(new_segs)
        new_segs = []
        for seg in temp:
            for v in H.neighbors(seg):
                if v not in source_segs and v not in K.nodes():
                    new_segs.append(v)
                if v in source_segs:
                    K_sources.add(v)
                K.add_edge(seg, v, weight=G[seg][v]['weight'])  
    for s in source_segs:
        if s in sink_segs:
            K.add_node(s)
            K_sources.add(s)
            K_sinks.add(s)
    if not K_sources:
        return None, None, None
    #print("SUCCESS. There are %d nanowires in the cluster out of all %d nanowires." % (len(K.nodes()), len(segs)))
    return K, K_sources, K_sinks