from Graph import *

# Method 1: check if the angle between current and previous vertex is too different from first angle in sequence
# Method 2: check if the angle between current and first vertex in the sequence is too different from first angle in sequence
# Method 3: method 2, but the max angle is affected by the distance from the first vertex in the sequence
# Method 4: method 1 and 2 are both active

def create_azimuth_graph(df, max_angle_change = 45, method = 3, 
                 max_angle_change_min = 5, distance_multiplier = None, splitting_method = 2):
    
    mean_distance = df["change"].mean()
    skip_next = False
    
    # Splitting distance is defined by a multiple of the mean
    if splitting_method == 1:
        if distance_multiplier == None: distance_multiplier = 5
        max_distance = mean_distance*distance_multiplier
    
    # Splitting distance is defined by the quantiles
    elif splitting_method == 2:
        if distance_multiplier == None: distance_multiplier = 1.5
        q25 = df["change"].quantile(0.25)
        q75 = df["change"].quantile(0.75)
        max_distance = (q75 + distance_multiplier*(q75 + q25))
        
    graph = Graph()
    
    # The start of the current sequence
    start_latitude = df.iloc[0, 0]
    start_longitude = df.iloc[0, 1]
    altitude = df.iloc[0, 2]
    
    # The last point
    start_vertex = Vertex(start_latitude, start_longitude, altitude, graph)
    last_latitude = df.iloc[1, 0]
    last_longitude = df.iloc[1, 1]
    last_altitude = df.iloc[1, 2]
    
    # The last angle to the vertex before it
    compare_last_angle = start_vertex.angle_to_point(last_latitude, last_longitude)
    
    # The last angle to the start vertex in the current_sequence
    compare_start_angle = start_vertex.angle_to_point(last_latitude, last_longitude)
    
    # For every point in the dataframe (in order)
    for i in range(2,len(df)):
        cur_latitude = df.iloc[i, 0]
        cur_longitude = df.iloc[i, 1]   
        cur_altitude = df.iloc[i, 2]
        
        last_distance = sqrt(pow(last_latitude-cur_latitude, 2) + pow(last_longitude-cur_longitude,2))

        # If the same coordinate is repeated
        if last_distance == 0:
           pass
        
        # If the last iteration noticed a skip, create a new start
        if skip_next == True:
            new_vertex = Vertex(last_latitude, last_longitude, last_altitude, graph)
            
            compare_last_angle = calculate_angle(last_latitude, cur_latitude, last_longitude, cur_longitude)
            compare_start_angle = compare_last_angle
            start_vertex = new_vertex

            start_latitude = last_latitude
            start_longitude = last_longitude
            
            skip_next = False            
        
        # If the distance is to large do not connect
        elif last_distance > max_distance:
            new_vertex = Vertex(last_latitude, last_longitude, last_altitude, graph)
            Edge([start_vertex, new_vertex])
            
            skip_next = True
        
        # Normal case of connecting
        else:
            if method == 1:
                last_angle = calculate_angle(last_latitude, cur_latitude, last_longitude, cur_longitude)

                condition = abs((abs(compare_last_angle-last_angle) + 180) % 360 - 180) > max_angle_change

            elif method == 2:
                start_angle = calculate_angle(start_latitude, cur_latitude, start_longitude, cur_longitude)
                
                condition = abs((abs(compare_start_angle-start_angle) + 180) % 360 - 180) > max_angle_change

            elif method == 3:
                start_angle = calculate_angle(start_latitude, cur_latitude, start_longitude, cur_longitude)
                
                start_distance = sqrt(pow(start_latitude-cur_latitude, 2) + pow(start_longitude-cur_longitude,2))
                max_angle = max(max_angle_change*mean_distance/start_distance, max_angle_change_min)
                
                condition = abs((abs(compare_start_angle-start_angle) + 180) % 360 - 180) > max_angle
            
            elif method == 4:
                last_angle = calculate_angle(last_latitude, cur_latitude, last_longitude, cur_longitude)
                start_angle = calculate_angle(start_latitude, cur_latitude, start_longitude, cur_longitude)
                
                condition = (abs((abs(compare_last_angle-last_angle) + 180) % 360 - 180) > max_angle_change)\
                or (abs((abs(compare_start_angle-start_angle) + 180) % 360 - 180) > max_angle_change)
            
            # Create new vertex if the angle is too different
            if condition:
                new_vertex = Vertex(last_latitude, last_longitude, last_altitude, graph)
                Edge([start_vertex, new_vertex])

                try: compare_last_angle = last_angle
                except:pass
                    
                try: compare_start_angle = start_angle
                except: pass

                start_vertex = new_vertex

                start_latitude = last_latitude
                start_longitude = last_longitude

        last_latitude = cur_latitude
        last_longitude = cur_longitude
        last_altitude = cur_altitude

    new_vertex = Vertex(last_latitude, last_longitude, last_altitude, graph)
    Edge([start_vertex, new_vertex])
    
    return graph