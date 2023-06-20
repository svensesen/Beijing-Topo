#Cell and Grid classes for dividing up the city into a grid
CITY_BORDERS = [40.19, 39.65, 115.98, 116.74]
from Graph import *

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
    
    def vertexize(self, threshold = 600):
        '''makes the Cell become a vertex if has more locations that threshold'''
        
        if self.number_of_locations > threshold:
            self.is_vertex = True

    

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
        '''for each cell number in the list, update the corresponding cell details'''
        
        for number in cell_numbers:
            self.cell_list[int(number)].number_of_locations += 1
            
    
    def vertexize(self, threshold):
        '''runs vertexize on each sell of the grid'''
        
        for cell in self.cell_list:
            cell.vertexize(threshold)
            
    def get_list_of_vertices(self) -> list:
        '''returns the numbers of all cells that are vertices (is_vertex() == True)'''
        
        return [i for i in range(len(self.cell_list)) if self.cell_list[i].is_vertex]
        
        
        
        
def latitude_to_row(latitude: float, latitude_cells_number: int, 
                    latitude_cell_length: float, CITY_BORDERS=CITY_BORDERS):
    '''For a given latitude in the row of a datafame, checks which row number in the grid the point belongs to'''
    
    for i in range(latitude_cells_number):
        if latitude < CITY_BORDERS[1] + i*latitude_cell_length: #we do not have to check the other containment!
            return i
    

def longitude_to_column(longitude: float, longitude_cells_number: int, 
                        longitude_cell_length: float, CITY_BORDERS=CITY_BORDERS):
    '''For a given longitude in the row of a datafame, checks which column number in the grid the point belongs to'''
    
    for i in range(longitude_cells_number):
        if longitude < CITY_BORDERS[2] + i*longitude_cell_length: #we do not have to check the other containment!
            return i
        

def row_column_to_cell_number(row_nr: int, column_nr: int, longitude_cells_nr: int):
    '''given the x and y in the grid, calculates the cell number'''
    
    return row_nr*longitude_cells_nr + column_nr
# lambda x,y,z: x*y+z


def find_subsequent_vertex_pairs(df_grouped, grid) -> dict:
    '''if a trajectory has visited vertices: A, A, A, B, C, D, D; we will obtain
    (A, B), (B, C), (C, D) as the key of the resulting dictionary'''

    dic = {}
    
    #For each trajectory, get the list of unique cells it traversed:
    for element in df_grouped['cell_number'].unique():
        
        #Proceed if there was more than 1 cell that was traversed:
        if len(element)>1:
            
            #If the person went from a vertex to a vertex and not just ANY cell:
            if grid.cell_list[int(element[0])].is_vertex and grid.cell_list[int(element[1])].is_vertex:
                
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
    
    #return g

    #Create all Edge objects - they will have the vertices added (and vertices will have edges added):
    for pair in edges:
        #v1 is the Vertex such that Vertex.ID = pair[0], same for v2
        try:
            v1, v2 = [v for v in g.vertices if v.id == pair[0]][0], [v for v in g.vertices if v.id == pair[1]][0]
        except:
            print("broken pair: ", pair)
        
        #Adding vertices to edges:
        e = Edge(vertices=set( (v1, v2) ) )
        
#         #Adding edges to vertices: #AT THIS POINT I ASSUME GRAPH GOT UPDATED BY THE EDGES TOO?
#         v1.add_edges(e)
#         v2.add_edges(e)
            
    return g