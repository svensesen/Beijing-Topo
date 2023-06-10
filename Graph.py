# Note each function starting with _ should NOT be manually invoked
from warnings import warn
from math import atan2, pi, sin, cos, sqrt, pow, degrees
from collections import deque
import matplotlib.pyplot as plt

import pickle

graph_counter = 0
vertex_counter = 0
edge_counter = 0

class Graph:
    def __init__(self, warnings = False):
        global graph_counter
        self.id = graph_counter
        graph_counter += 1

        self.warnings = warnings

        self.vertices = set()
        self.edges = set()

    def __repr__(self):
        return "<Graph>"
    
    def __hash__(self):
        return self.id
    
    def combine_graph(self, graph):
        self.vertices.update(graph.vertices)
        self.edges.update(graph.edges)

    def merge_close_vertices(self, max_distance = 0.0001):
        queue = deque(self.vertices)

        # For each vertex we check
        while len(queue) != 0:
            cur_vertex = queue.popleft()
            to_combine = set()

            # For each other vertex check if it is close
            for other_vertex in self.vertices:
                if other_vertex != cur_vertex:
                    distance = cur_vertex.short_distance_to_vertex(other_vertex)
                    if distance <= max_distance:
                        to_combine.add(other_vertex)
            
            # If there were close vertices combine them
            if len(to_combine) > 0:
                latitudes = [vertex.latitude for vertex in to_combine] + [vertex.latitude]
                longitudes = [vertex.longitude for vertex in to_combine] + [vertex.longitude]

                cur_vertex.latitude = sum(latitudes)/len(latitudes)
                cur_vertex.longitude = sum(longitudes)/len(longitudes)

                for vertex in to_combine: 
                    for edge in vertex.edges:
                        edge.remove_vertices(vertex)
                        edge.add_vertices(cur_vertex)
                    queue.remove(vertex)

                self.vertices = self.vertices.difference(to_combine)
    
    def merge_vertices_to_edges(self):
        for cur_vertex in self.vertices:
            if 


    
    def add_vertices(self, vertices, called = False):
        if not isinstance(vertices, set):
            try:
                vertices = set(vertices)
            except:
                vertices = {vertices}
        
        for value in vertices:
            if not isinstance(value, Vertex):
                raise TypeError(f"Input for vertices contained a non-vertex, namely {type(value)}")
        
        if self.warnings and self.vertices.intersection(vertices):
            warn("Attempting to add vertices to graph that are already there")

        self.vertices.update(vertices)

        for vertex in vertices:
            for edge in vertex.edges:
                self.edges.add(edge)
        
        if not called:
            vertex.set_graph(self, called = True)
    
    def remove_vertices(self, vertices, called = False):
        if not isinstance(vertices, set):
            try:
                vertices = set(vertices)
            except:
                vertices = {vertices}

        for value in vertices:
            if not isinstance(value, Vertex):
                raise TypeError(f"Input for vertices contained a non-vertex, namely {type(value)}")

        if self.warnings and not vertices.issubset(self.vertices):
            warn("Attempting to remove vertices from the graph that are not there")

        self.vertices = self.vertices.difference(vertices)

        for vertex in vertices:
            for edge in vertex.edges:
                if not self.vertices.intersection(edge.vertices):
                    self.edges.remove(edge)
        
        if not called:
            vertex.set_graph(None, called = True)

    def closest_vertex_to_point(self, latitude, longitude):        
        closest_vertex = None
        closest_distance = float("inf")
        for vertex in self.graph.vertices:
            distance = vertex.distance_to_point(latitude, longitude)

            if distance < closest_distance:
                closest_vertex = vertex
                closest_distance = distance

        return closest_vertex
    
    def vertices_in_square(self, latitude0, latitude1, longitude0, longitude1):
        if latitude1 < latitude0:
            latitude0, latitude1 = latitude1, latitude0

        if longitude1 < longitude0:
            longitude0, longitude1 = longitude1, longitude0

        found_vertices = set()
        for vertex in self.graph.vertices:
            if (vertex.latitude >= latitude0) and (vertex.latitude <= latitude1) and \
            (vertex.longitude >= longitude0) and (vertex.longitude <= longitude1):
                found_vertices.add(vertex)
        
        return found_vertices

    def _add_edges(self, edges):
        self.edges.update(edges)
    
    def _remove_edges(self, edges):
        self.edges = self.edges.difference(edges)
        
    def remove_unconnected_vertices(self):
        self.vertices = {vertex for vertex in self.vertices if len(vertex.edges) != 0}
        
    def save(self, filename: str): #works, but not loading yet
        with open(filename, 'wb') as outp:
            pickle.dump(self, outp, -1)

    
class Vertex:
    def __init__(self, latitude, longitude, altitude, graph = None, edges = None, warnings = False):
        global vertex_counter
        self.id = vertex_counter
        vertex_counter += 1

        self.warnings = warnings

        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude

        self.graph = None

        self.edges = set()
        if edges != None:
            self.set_edges(edges)

        if graph != None:
            self.set_graph(graph)

    def __repr__(self): 
        return f"<Vertex: {self.latitude},{self.longitude},{self.altitude}>"
    
    def __hash__(self):
        return self.id
    
    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return (self.latitude + self.longitude) < (other.latitude + other.longitude)

    def __le__(self, other):
        return (self.latitude + self.longitude) <= (other.latitude + other.longitude)
    
    def __gt__(self, other):
        return (self.latitude + self.longitude) > (other.latitude + other.longitude)

    def __ge__(self, other):
        return (self.latitude + self.longitude) >= (other.latitude + other.longitude)
    
    def add_edges(self, edges, called = False):
        if not isinstance(edges, set):
            try:
                edges = set(edges)
            except:
                edges = {edges}
    
        for value in edges:
            if not isinstance(value, Edge):
                raise TypeError(f"Input for edges contained a non-Edge, namely {type(value)}")
        
        if self.warnings and self.edges.intersection(edges):
            warn("Attempting to add edge(s) to vertex that are already there")

        self.edges.update(edges)

        if self.graph != None:
            self.graph._add_edges(edges)
        
        if not called:
            for edge in edges:
                edge.add_vertices(self, called = True)
        
        for edge in self.edges:
            edge._set_graph(self.graph)
    
    def remove_edges(self, edges, called = False):
        if not isinstance(edges, set):
            try:
                edges = set(edges)
            except:
                edges = {edges}
    
        for value in edges:
            if not isinstance(value, Edge):
                raise TypeError(f"Input for edges contained a non-Edge, namely {type(value)}")
        
        if self.warnings and not edges.issubset(self.edges):
            warn("Attempting to remove edges(s) from the vertex that are not there")

        self.edges = self.edges.difference(edges)

        if self.graph != None:
            self.graph._remove_edges(edges)
        
        if not called:
            for edge in edges:
                edge.remove_vertices(self, called = True)
    
        for edge in self.edges:
            if not edge.vertices:
                edge._set_graph(None)
    
    def set_graph(self, graph, called = False):
        if (not isinstance(graph, Graph)) and (graph != None):
            raise TypeError(f"A graph must be a Graph or None, got {type(graph)}")
        
        if not called:
            if self.graph != None:
                self.graph.remove_vertices(self, called = True)
            
            if graph != None:
                graph.add_vertices(self, called = True)
        
        self.graph = graph

        for edge in self.edges:
            edge._set_graph(self.graph)

    def angle_to_vertex(self, vertex):
        return calculate_angle(self.latitude, vertex.latitude, self.longitude, vertex.longitude)
    
    def angle_to_point(self, latitude, longitude):
        return calculate_angle(self.latitude, latitude, self.longitude, longitude)
    
    def distance_to_vertex(self, vertex):
        return calculate_distance(self.latitude, vertex.latitude, self.longitude, vertex.longitude)
    
    def distance_to_point(self, latitude, longitude):
        return calculate_distance(self.latitude, latitude, self.longitude, longitude)
    
    def short_angle_to_vertex(self, vertex):
        return short_angle(self.latitude, vertex.latitude, self.longitude, vertex.longitude)
    
    def short_angle_to_point(self, latitude, longitude):
        return short_angle(self.latitude, latitude, self.longitude, longitude)
    
    def short_distance_to_vertex(self, vertex):
        return short_distance(self.latitude, vertex.latitude, self.longitude, vertex.longitude)
    
    def short_to_point(self, latitude, longitude):
        return short_distance(self.latitude, latitude, self.longitude, longitude)
    
    def closest_vertex(self):
        if self.graph == None:
            return None
        
        closest_vertex = None
        closest_distance = float("inf")
        for vertex in self.graph.vertices:
            distance = self.distance_to_vertex(vertex)

            if distance < closest_distance:
                closest_vertex = vertex
                closest_distance = distance

        return closest_vertex 

    def short_distance_to_edge(self, edge):
        if len(edge.vertices) != 2:
            return -1
        
        v0 = min(edge.vertices)
        v1 = max(edge.vertices)
        return calculate_distance_point_line(v0.latitude, v1.latitude, v0.longitude, v1.longitude,\
                                             self.latitude, self.longitude)

class Edge:
    def __init__(self, vertices = None, warnings = False):
        global edge_counter
        self.id = edge_counter
        edge_counter += 1

        self.warnings = warnings

        self.vertices = set()
        self.graph = None
        
        if vertices != None:
            self.add_vertices(vertices)
        
    def __repr__(self):
        if len(self.vertices) == 0:
            return "<Edge>"
        
        elif len(self.vertices) == 1:
            return f"<Edge: ({next(iter(self.vertices))}) "
        
        elif len(self.vertices) == 2:
            return f"<Edge: ({min(self.vertices)}) -> ({max(self.vertices)})>"
    
    def __hash__(self):
        return self.id
        
    def add_vertices(self, vertices, called = False):
        if not isinstance(vertices, set):
            try:
                vertices = set(vertices)
            except:
                vertices = {vertices}

        for value in vertices:
            if not isinstance(value, Vertex):
                raise TypeError(f"Input for vertices contained a non-Vertex, namely {type(value)}")
        
        if len(self.vertices) + len(vertices) > 2:
            raise ValueError(f"An edge may have at most 2 vertices, has {len(self.vertices)}, adding {len(vertices)}")
        
        if self.warnings and self.vertices.intersection(vertices):
            warn("Attempting to add vertices to edge that are already there")

        self.vertices.update(vertices)
        
        if not called:
            for vertex in vertices:
                vertex.add_edges(self, called = True)
    
    def remove_vertices(self, vertices, called = False):
        if not isinstance(vertices, set):
            try:
                vertices = set(vertices)
            except:
                vertices = {vertices}

        for value in vertices:
            if not isinstance(value, Vertex):
                raise TypeError(f"Input for vertices contained a non-Vertex, namely {type(value)}")
        
        if len(vertices) != 2:
            raise ValueError(f"An edge may have at most 2 vertices, has {len(self.vertices)}, adding {len(vertices)}")
        
        if self.warnings and not vertices.issubset(self.vertices):
            warn("Attempting to remove vertices from the edge that are not there")

        self.vertices = self.vertices.difference(vertices)

        if not called:
            for vertex in vertices:
                vertex.remove_edges(set(self), called = True)
    
    def _set_graph(self, graph):
        self.graph = graph

    def short_distance_to_vertex(self, vertex):
        if len(self.vertices) != 2:
            return -1

        v0 = min(self.vertices)
        v1 = max(self.vertices)
        return calculate_distance_point_line(v0.latitude, v1.latitude, v0.longitude, v1.longitude,\
                                             vertex.latitude, vertex.longitude)

    @property
    def length(self):
        if len(self.vertices) != 2:
            return -1

        v0 = min(self.vertices)
        v1 = max(self.vertices)

        return calculate_distance(v0.latitude, v1.latitude, v0.longitude, v1.longitude)
        
    @property
    def angle(self):
        if len(self.vertices) != 2:
            return -1
        
        v0 = min(self.vertices)
        v1 = max(self.vertices)

        return calculate_angle(v0.latitude, v1.latitude, v0.longitude, v1.longitude)


def calculate_distance(latitude0, latitude1, longitude0, longitude1):
        R = 6371009
        ph0 = latitude0 * pi/180
        ph1 = latitude1 * pi/180
        dlat = (latitude1 - latitude0) * pi/180
        dlong = (longitude1 - longitude0) * pi/180

        a = sin(dlat/2) * sin(dlat/2) + cos(ph0) * cos(ph1) * sin(dlong/2) * sin(dlong/2)
        c = 2 * atan2(sqrt(a), sqrt(1-a))

        d = R * c

        return d


def calculate_angle(latitude0, latitude1, longitude0, longitude1):
    y = sin(longitude1-longitude0) * cos(latitude1)
    x = cos(latitude0) * sin(latitude1) - sin(latitude0) * cos(latitude1) * cos(longitude1 - longitude0)
    ph = atan2(y, x)
    b = (ph*180/pi + 360) % 360

    return b


def short_angle(latitude0, latitude1, longitude0, longitude1):
    y = abs(latitude0 - latitude1)
    x = abs(longitude0 - longitude1)
    b = degrees(atan2(y, x))

    return b

def short_distance(x0, x1, y0, y1):
    return sqrt(pow(x1-x0, 2) + pow(y1-y0, 2))
    

def calculate_distance_point_line(x0, x1, y0, y1, px, py):
    return abs((x1-x0)*(y0-py) - (x0-px)*(y1-y0))/sqrt(pow((x1-x0), 2) + pow((y1-y0), 2))


def plot_sus(lis, base = None, marker = ',', alti_group: bool = False):
    #lis = a list of Vertex and Edge objects
    #base = the base image to plot over. path to a png (preferably)
    #marker = the marker to be used while plotting, default = pixel
    #alti_group = whether you wanna color by altitude or not (bool)
    get_line_x = lambda edge: np.array([(list(edge.vertices)[i][0], list(edge.vertices)[i+1][0]) for i in range(0, len(edge.vertices)-1)])
    get_line_y = lambda edge: np.array([(list(edge.vertices)[i][1], list(edge.vertices)[i+1][1]) for i in range(0, len(edge.vertices)-1)])

    if base:
       im = plt.imread(base) #needs to be path to image, preferably png
       implot = plt.imshow(im) 
    if alti_group:
        dict_colors = {}
        alt_groups = set([i.altitude for i in vertex_list])
        cntr = 0
    
        for group in alt_groups:
            dict_colors[group] = mcolors.cnames.keys()[cntr]
            cntr += 1
    lat, lon, alt  = [], [], []
    for v in lis:
        if type(v) == Vertex:
           lat.append(v.latitude)
           lon.append(v.longitude)
           alt.append(v.altitude)
       #if alt:
           #plt.scatter(lat, lon, alt, marker=marker)
       #else:
           if alti_group:
               plt.scatter(lat, lon, c = alt.map(dict_colors), marker = marker)
               plt.legend()
           else:
               plt.scatter(lat, lon, marker = marker)
        elif type(v) == Edge:
            x = get_line_x(v)
            y = get_line_y(v)
            
            for i in range(len(x)):
                plt.plot(x[i], y[i], 'ro-')