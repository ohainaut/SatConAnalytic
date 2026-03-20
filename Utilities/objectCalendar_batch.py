#!/usr/bin/env python3


import numpy as np
import objectCalendar

import objectCalendar

myStars = np.array([ [17., 20.], [18., 10],
    [19, -10], [20,12], [21, -43],[22, -4],[23, -64],[0,-14],[1,-32],[2,-13.],
    [20.44,-14], [22.1,-20],[1,4],[2,20],
    [23.1,-60],
    [21.39,-11.], [2., -61.]
    ])

myStars = myStars* np.array([-15.,1.])  # convert hours to deg
myStars = myStars+ np.array([90.,0.])  # shift to night side


for star in myStars:
    ra = star[0] + np.random.uniform(-0.5,0.5)
    dec = star[1] + np.random.uniform(-0.5,0.5)
    print(f'Star RA={ra:.2f} deg Dec={dec:.2f} deg')


    # Define arguments as you would on command line
    args = [
        '-a', f'{ra}', 
        '-d', f'{dec}',
        '-C', 'SLOWGWAK',#'TODAY', #'ALL', #'SL1old',        # Starlink constellation
        '-T', 'FORSimg',
        '--mode', 'EFFECT'
    ]

    print(args)

    objectCalendar.main(args)
    print("Single star simulation completed")

print("All simulations completed!")
