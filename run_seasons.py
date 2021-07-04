import numpy as np
from msgp import MSGP_seasons
from utils import print_results_seasons


data = np.load('AirQualityUCI.npy')
# Keep only 6 first attributes
data = data.reshape((389, 24, 13))[:, :, :6]

minSup = int(389*0.15)
print_results_seasons(MSGP_seasons(data, minSup), 24, 6)
