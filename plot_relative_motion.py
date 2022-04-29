from posixpath import split
from random import sample
import pygmt
import os
import pandas as pd
import numpy as np

def relative_motion(lambda_x, lambda_p, phi_x, phi_p, omega):
    #########################################################
    # We can use this function to calculate the relative 
    # motion of the specific plate.
    # Please note the unit of each parameter!
    # lambda_x, lambda_p, phi_x, phi_p are in degree
    # omega -> [degree/year]
    # It refer to the example from textbook p.21
    #########################################################
    radius = 6371. * 1e+05              # From [km] to [cm] 
    lambda_x = np.deg2rad(lambda_x)
    lambda_p = np.deg2rad(lambda_p)
    phi_x = np.deg2rad(phi_x)
    phi_p = np.deg2rad(phi_p)
    omega = np.deg2rad(omega)

    angle_a = np.arccos(np.sin(lambda_x)*np.sin(lambda_p)\
        +np.cos(lambda_x)*np.cos(lambda_p)*np.cos(phi_p-phi_x))
    
    angle_C = np.arcsin(np.cos(lambda_p) * np.sin(phi_p-phi_x) / np.sin(angle_a))
    velocity = omega * radius * np.sin(angle_a)
    angle_C = np.rad2deg(angle_C)
    beta = 90 + angle_C
    return angle_C, velocity, beta

def read_the_boundary_file(file_for_reading):
    f = open(file_for_reading, 'r')
    data = f.readlines()
    f.close()
    return data

def find_the_header(data):
    header_index = []
    for ii, line in enumerate(data):
        if line[:3] == '***':
            header_index.append(ii)
    return header_index

def get_data_for_each_plate_boundary(data, line1, line2):
    lon_list, lat_list = [], []
    for line in data[line1:line2]:
        line = line.strip()
        info = line.split(',')
        lon, lat = float(info[0]), float(info[1])
        lon_list.append(lon)
        lat_list.append(lat)
    return lon_list, lat_list


# Download the grid file (for topography)
region = [118, 127, 20, 26]
grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)

# Create your basemap and plot the coastal line
fig = pygmt.Figure()
fig.grdimage(grid, region=region, projection='M6i', cmap='ETOPO1.cpt')
fig.coast(region=region, projection="M6i",
            shorelines=['0.2p', 'black'], frame=['a2f0.5', 'WSne'], area_thresh=0.1)


# Read the file with the boundary data near Taiwan
data = read_the_boundary_file('PB2002_boundaries.dig.txt')
header_index = find_the_header(data)
for ii, segment in enumerate(header_index):
    if ii == 0:
        line1 = 1
        line2 = segment
    else:
        line1 = line2+2
        line2 = segment
    lon_list, lat_list = get_data_for_each_plate_boundary(data, line1, line2)
    boundary_df = pd.DataFrame(zip(lon_list, lat_list), columns=['lon', 'lat'])
    fig.plot(x=boundary_df.lon, y=boundary_df.lat, pen="0.8p,red") 


### Create the sample points 
### (It's not necessary in your homework because the question has specify the points you need to calculate)
sample_points_lon, sample_points_lat = [122, 123, 124], [22, 22, 22]
sample_points_df = pd.DataFrame(zip(sample_points_lon, sample_points_lat), columns=['lon', 'lat'])
fig.plot(x=sample_points_df.lon, y=sample_points_df.lat, style="c0.1c", color='red', pen="black")

### Use the function "relative_motion" to calculate, please enter your points for calculating and the other parameters
angle_C, velocity, beta = relative_motion(lambda_x=sample_points_df.lon,
                                          lambda_p=158.9, 
                                          phi_x=sample_points_df.lat, 
                                          phi_p=46.67, 
                                          omega=1.035e-06)
beta = 360 - beta ## If the angle of plate A <--> plate B is beta, and the angle of Plate B <--> Plate A is (360-beta)

### make the velocity project into 2 component N and E 
vec_N = velocity * np.cos(np.deg2rad(beta))
vec_E = velocity * np.sin(np.deg2rad(beta))

### plot the velocity arrows
points_num = len(sample_points_df)
df = pd.DataFrame(
    data={
        "x": sample_points_df.lon,
        "y": sample_points_df.lat,
        "east_velocity": vec_E,
        "north_velocity": vec_N,
        "east_sigma": np.zeros(points_num),
        "north_sigma": np.zeros(points_num),
        "correlation_EN": np.zeros(points_num)
    }
)
fig.velo(
    data=df,
    pen="1p,black",
    spec="e0.2/0",
    color='red'
)

#################################################################################
# Here we plot the reference arrow with 7.5 cm/yr
df_ex = pd.DataFrame(
    data={
        "x": [120],
        "y": [21.1],
        "east_velocity": [-7.5],
        "north_velocity": [0],
        "east_sigma": [0],
        "north_sigma": [0],
        "correlation_EN": [0]
    }
)
fig.velo(
    data=df_ex,
    pen="1p,black",
    spec="e0.2/0",
    color='red'
)
# Plot the text as label
fig.text(x=124.5, y=21.5, text='Philippine Sea Plate', font='18p,Helvetica-Bold,black')
fig.text(x=119.5, y=25, text='Eurasian Plate', font='18p,Helvetica-Bold,black')
fig.text(x=119.5, y=20.85, text='7.5 cm/yr', font='10p')
#################################################################################
fig.show()
fig.savefig('Taiwan.png')