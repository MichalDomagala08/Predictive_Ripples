import numpy as np

screen_width = 1920
mean_viewing_distance = 590
screen_width_mm = 530
deg2_px_radius = mean_viewing_distance * np.tan(np.radians(2)) * screen_width / screen_width_mm
print(deg2_px_radius)