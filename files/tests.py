import numpy as np
import matplotlib.pyplot as plt


def read_from_file(filename):
    file = np.genfromtxt(filename, skip_header=18, skip_footer=1)
    x = file[:1]
    y = file[:2]
    return (x,y)

read_from_file("laser_he_ne_0.txt")