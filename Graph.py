# Note each function starting with _ should NOT be manually invoked
from warnings import warn

graph_counter = 0
vertex_counter = 0
edge_counter = 0

class Graph:
    def __init__(self, warnings = False):
        self.vertices = set()
        self.edges = set()
        self.warnings = warnings

        global graph_counter
        self.id = graph_counter
        graph_counter += 1

    def __repr__(self):
        return "<Graph>"
    
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

    def _add_edges(self, edges):
        self.edges.update(edges)
    
    def _remove_edges(self, edges):
        self.edges = self.edges.difference(edges)


class Vertex:
    def __init__(self, longitude, latitude, altitude, graph = None, edges = None, warnings = False):
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude

        self.graph = graph

        self.warnings = warnings

        self.edges = set()
        if edges != None:
            self.set_edges(edges)
        
        global vertex_counter
        self.id = vertex_counter
        vertex_counter += 1
    
    def __repr__(self): 
        return f"<Vertex: {self.longitude},{self.latitude},{self.altitude}>"
    
    def __hash__(self):
        return self.id
    
    def __eq__(self, other):
        return (self.longitude == other.longitude) and (self.latitude == other.latitude)

    def __lt__(self, other):
        return (self.longitude + self.latitude) < (other.longitude + other.latitude)

    def __le__(self, other):
        return (self.longitude + self.latitude) <= (other.longitude + other.latitude)
    
    def __gt__(self, other):
        return (self.longitude + self.latitude) > (other.longitude + other.latitude)

    def __ge__(self, other):
        return (self.longitude + self.latitude) >= (other.longitude + other.latitude)
    
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
                self.graph.add_vertices(self, called = True)
        
        self.graph = graph

        for edge in self.edges:
            edge._set_graph(self.graph)

        

class Edge:
    def __init__(self, graph = None, vertices = None, warnings = False):
        global edge_counter
        self.id = edge_counter
        edge_counter += 1

        self.vertices = set()
        self.graph = graph
        self.warnings = warnings

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
        
        if len(self.vertices) + len(vertices) != 2:
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
    

        
