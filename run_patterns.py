import numpy as np
from msgp import MSGP_patterns
from utils import print_results_patterns


data = np.load('AirQualityUCI.npy')
# Keep only 6 first attributes
data = data.reshape((389, 24, 13))[:, :, :6].reshape((389*24, 6))

minSup = int(389*0.15)
print_results_patterns(MSGP_patterns(data, 24, minSup), 24, 6)
