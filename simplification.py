from Graph import *

def simplify_double_edges(g: Graph):
    '''if there are 2 edges within a pair of vertices, leave only 1
    This is possible due to perhaps multiple edge objects (of different IDs) connecting the same Vertices
    
    Also fixes broken edges'''
    
    unique_vertex_connections = []
    edges_ids_to_be_deleted = []

    for edge in g.edges:
        
        #Not broken edge (2 endpoints)
        pair = (list(edge.vertices)[0].id, list(edge.vertices)[1].id)

        if pair not in unique_vertex_connections:
            unique_vertex_connections.append(pair)
        else:
            edges_ids_to_be_deleted.append(edge.id)
    
    print(len(edges_ids_to_be_deleted))
    print(edges_ids_to_be_deleted)
    g.edges = set([edge for edge in g.edges if edge.id not in edges_ids_to_be_deleted])
    
    for vertex in g.vertices:
        vertex.edges = set([edge for edge in vertex.edges if edge.id not in edges_ids_to_be_deleted])
    
    return g
	
	
	
def delete_chain_vertices(g: Graph, threshold = 0.05):
    '''The threshold is not very thought-through btw'''
    
    # 1) for every vertex:
    for vertex in g.vertices:
    
    # 2) if vertex has precisely 2 edges:
        if len(vertex.edges) == 2:
            edges = list(vertex.edges)
            
            # 3) if their angles are matching (below some threshold):
            if abs(list(vertex.edges)[0].angle % 180 - list(vertex.edges)[1].angle % 180) < threshold:
                
                # 4) Find all vertices this chain connects
                #That means, we'll get the outer vertices twice and the middle vertex once
                all_vertices = list(edges[0].vertices) + list(edges[1].vertices)
                
                # 5) Sometime, due to completely fucked up errors, we will still encounter double edges instead of chains
                # This is why we use the Counter (histogram) to check for that cases
                counter = Counter(all_vertices)
                #print(counter)
                
                # 6) To refer to the outer vertices, check which appear in the Counter once (middle will be twice)
                outer_vertices = [out for out in counter if counter[out] == 1]
                
                # 7) If there are no outer vertices, we have found a redundant edge. Hence, remove any edge of the two:
                if len(outer_vertices) == 0:
                    g.edges.remove(edges[0])
                
                # 8) In case this is indeed a chain:
                else:
                    
                    # 9) join the other 2 vertices using an edge
                    new_edge = {Edge(vertices=set(outer_vertices))}

                    # 10) remove the middle vertex
                    g.remove_vertices(vertex)

    return g