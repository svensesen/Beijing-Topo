#Cell and Grid classes for dividing up the city into a grid
CITY_BORDERS = [40.19, 39.65, 115.98, 116.74]
from Graph import *
from collections import Counter

class Cell():
    '''Just your regular Square kind of a class'''
    def __init__(self, top: float, left: float, bottom: float, right: float):
        self.top = top
        self.left = left
        self.bottom = bottom
        self.right = right
        self.number_of_locations = 0
        self.is_vertex = False
        
    def __str__(self):
        return "Cell: top = " + str(round(self.top, 5)) + ", left = " + str(round(self.left, 5)) + \
    ", bottom = " + str(round(self.bottom, 5)) + ", right = " + str(round(self.right, 5))

    def belongs(self, latitude: float, longitude: float):
        '''checks if a point belongs to the cell - pass the point as (y, x)'''
        return (latitude > self.bottom and latitude <= self.top and longitude < self.right and longitude >= self.left)

    

class Grid():
    '''pretty much just a list of Cells that is available at self.cell_list'''
    
    def __init__(self, latitude_cells_nr: int, longitude_cells_nr: int, city_borders = CITY_BORDERS):
        self.latitude_cells_nr = latitude_cells_nr
        self.longitude_cells_nr = longitude_cells_nr
        self.city_borders = city_borders
        
        self.latitude_cell_length = round((city_borders[0] - city_borders[1])/latitude_cells_nr, 5)
        self.longitude_cell_length = round((city_borders[3] - city_borders[2])/longitude_cells_nr, 5)
        
        self.cell_list = []
        for i in range(latitude_cells_nr):
            for j in range(longitude_cells_nr):
                self.cell_list.append(Cell(top=city_borders[1] + self.latitude_cell_length*i,
                                     left=city_borders[2] + self.longitude_cell_length*j,
                                     bottom=city_borders[1] + self.latitude_cell_length*(i+1),
                                     right=city_borders[2] + self.longitude_cell_length*(j+1)))
                
    def __str__(self):
        return "Grid of side length (latitude) " + str(self.latitude_cell_length) + \
    " and (longitude) " + str(self.longitude_cell_length) + " and total number of cells: " + \
    str(self.latitude_cells_nr*self.longitude_cells_nr)

        
    def feed_list_of_points(self, points: list):
        '''point is a tuple (y, x) - DO NOT MESS THE Y AND X UP!
        increase the number_of_locations for each cell, if the point belongs there'''
        
        for point in points:
            for cell in self.cell_list:
                if cell.belongs(point[0], point[1]):
                    cell.number_of_locations += 1
                    break
                    
                    
    def feed_list_of_cell_numbers(self, cell_numbers: list):
        '''for each cell number in the list, update the corresponding cell details
        cell_numbers can also be a Pandas Series'''
        
        # A histogram of how many times was each cell number present:
        c = Counter(cell_numbers)
        
        # For each cell number, increase its number of locations by as much as the histogram value:
        for key in c:
            self.cell_list[key].number_of_locations += c[key]
            
    
    def vertexize(self, threshold):
        '''Compares the number of GPS locations in each cell with the threshold and flips the is_vertex flags is needed'''
        
        for cell in self.cell_list:
            if cell.number_of_locations > threshold:
                cell.is_vertex = True
            
    def get_list_of_vertices(self) -> list:
        '''returns the numbers of all cells that are vertices (is_vertex() == True)'''
        
        return [i for i in range(len(self.cell_list)) if self.cell_list[i].is_vertex]
        
        
        
def latitude_to_cell_nr(latitude: float, longitude_cells_nr: int, 
                    latitude_cell_length: float, CITY_BORDERS=CITY_BORDERS):
    '''For a given latitude in the row of a datafame, checks which row number in the grid the point belongs to'''
    
    #position - offset, divided by cell length, all times the number of cells in a row (not column!)
    return int(((latitude-CITY_BORDERS[1])//latitude_cell_length)*longitude_cells_nr)
    

def longitude_to_cell_nr(longitude: float, longitude_cell_length: float, CITY_BORDERS=CITY_BORDERS):
    '''For a given longitude in the row of a datafame, checks which column number in the grid the point belongs to'''
    
    #position - offset, divided by cell length
    return int((longitude-CITY_BORDERS[2])//longitude_cell_length)
        


def find_subsequent_vertex_pairs(df_grouped, grid) -> dict:
    '''if a trajectory has visited vertices: A, A, A, B, C, D, D; we will obtain
    (A, B), (B, C), (C, D) as the key of the resulting dictionary'''

    dic = {}
    
    #For each trajectory, get the list of unique cells it traversed:
    for element in df_grouped['cell_number'].unique():
        
        #Proceed if there was more than 1 cell that was traversed:
        if len(element)>1:
            
            #If the person went from a vertex to a vertex and not just ANY cell:
            if grid.cell_list[int(element[0])].is_vertex and grid.cell_list[int(element[-1])].is_vertex:
                
                #Get all pairs of subsequent vertices
                list_of_pairs = [(int(element[i]), int(element[i+1])) for i in range(len(element)-1)]
                
                #Histogram code:
                for pair in list_of_pairs:
                    if pair not in dic:
                        dic[pair] = 1
                    else:
                        dic[pair] += 1

    return dic
    
    
def rav_graph_to_sven_graph_2(grid, edges):
    ''''''
    
    #Create a Graph object:
    g = Graph()
    
    #Create all Vertex objects and add them to the Graph:
    vertex_ids = grid.get_list_of_vertices()
    for ID in vertex_ids:
        #Feel free to change the top-left into bottom-right or the middle coordinates!
        v = Vertex(latitude=grid.cell_list[ID].top, longitude=grid.cell_list[ID].left, altitude=100)
        v.id = ID
        g.add_vertices(v)
        
#     vertex_list = [Vertex(latitude=grid.cell_list[ID].top, longitude=grid.cell_list[ID].left, altitude=100) for ID in vertex_ids]
#     list(map(lambda x: x.id=, vertex_list))

    # The following is the fastest way to solve this, but for some reason returns some errors:
#     g.add_vertices([Vertex(latitude=grid.cell_list[vertex_id].top, longitude=grid.cell_list[vertex_id].left, \
#                            altitude=100, ID=vertex_id) for vertex_id in vertex_ids])

    #Create all Edge objects - they will have the vertices added (and vertices will have edges added):
    for pair in edges:
        #v1 is the Vertex such that Vertex.ID = pair[0], same for v2
        try:
            v1, v2 = [v for v in g.vertices if v.id == pair[0]][0], [v for v in g.vertices if v.id == pair[1]][0]
            #Adding vertices to edges:
            e = Edge(vertices=set( (v1, v2) ) )
        except:
            print("broken pair: ", pair)
        
    #AT THIS POINT I ASSUME GRAPH GOT UPDATED BY THE EDGES TOO
    return g