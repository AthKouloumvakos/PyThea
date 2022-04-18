'''
Configuration for the sliders
'''

############################################################
# Geometrical Models

sliders_Spheroid = {
    'height': {'Standard': {'min': 1., 'max': 30., 'step': 0.01, 'def': 2.},
               '<10Rsun': {'min': 1., 'max': 10., 'step': 0.01, 'def': 2.},
               '>10Rsun': {'min': 1., 'max': 30., 'step': 0.01, 'def': 2.},
               '>30Rsun': {'min': 1., 'max': 215., 'step': 0.1, 'def': 2.},
               'unit': '[R_sun]'},
    'kappa': {'Standard': {'min': 0., 'max': 2., 'step': 0.01, 'def': 1.},
              '<10Rsun': {'min': 0., 'max': 2., 'step': 0.01, 'def': 1.},
              '>10Rsun': {'min': 0., 'max': 2., 'step': 0.01, 'def': 1.},
              '>30Rsun': {'min': 0., 'max': 2., 'step': 0.01, 'def': 1.},
              'unit': ''},
    'epsilon': {'Standard': {'min': -1., 'max': 1., 'step': 0.01, 'def': 0.},
                '<10Rsun': {'min': -1., 'max': 1., 'step': 0.01, 'def': 0.},
                '>10Rsun': {'min': -1., 'max': 1., 'step': 0.01, 'def': 0.},
                '>30Rsun': {'min': -1., 'max': 1., 'step': 0.01, 'def': 0.},
                'unit': ''},
    'rcenter': {'Standard': {'min': 1., 'max': 30., 'step': 0.01, 'def': 1.},
                '<10Rsun': {'min': 1., 'max': 10., 'step': 0.01, 'def': 1.},
                '>10Rsun': {'min': 1., 'max': 30., 'step': 0.01, 'def': 1.},
                '>30Rsun': {'min': 1., 'max': 215., 'step': 0.1, 'def': 1.},
                'unit': '[R_sun]'},
    'radaxis': {'Standard': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                '<10Rsun': {'min': 0., 'max': 10., 'step': 0.01, 'def': 1.},
                '>10Rsun': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                '>30Rsun': {'min': 0., 'max': 215., 'step': 0.1, 'def': 1.},
                'unit': '[R_sun]'},
    'orthoaxis1': {'Standard': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                   '<10Rsun': {'min': 0., 'max': 10., 'step': 0.01, 'def': 1.},
                   '>10Rsun': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                   '>30Rsun': {'min': 0., 'max': 215., 'step': 0.1, 'def': 1.},
                   'unit': '[R_sun]'},
}

sliders_Ellipsoid = {
    'alpha': {'Standard': {'min': 0.5, 'max': 1.5, 'step': 0.01, 'def': 1.},
              '<10Rsun': {'min': 0.5, 'max': 1.5, 'step': 0.01, 'def': 1.},
              '>10Rsun': {'min': 0.5, 'max': 1.5, 'step': 0.01, 'def': 1.},
              '>30Rsun': {'min': 0.5, 'max': 1.5, 'step': 0.01, 'def': 1.},
              'unit': ''},
    'tilt': {'Standard': {'min': -90., 'max': 90., 'step': 0.25, 'def': 0.},
             '<10Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             '>10Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             '>30Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             'unit': '[deg.]'},
    'orthoaxis2': {'Standard': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                   '<10Rsun': {'min': 0., 'max': 10., 'step': 0.01, 'def': 1.},
                   '>10Rsun': {'min': 0., 'max': 30., 'step': 0.01, 'def': 1.},
                   '>30Rsun': {'min': 0., 'max': 215., 'step': 0.1, 'def': 1.},
                   'unit': '[R_sun]'}
}

sliders_GCS = {
    'height': {'Standard': {'min': 1., 'max': 30., 'step': 0.01, 'def': 4.},
               '<10Rsun': {'min': 1., 'max': 10., 'step': 0.01, 'def': 4.},
               '>10Rsun': {'min': 1., 'max': 30., 'step': 0.01, 'def': 4.},
               '>30Rsun': {'min': 1., 'max': 215., 'step': 0.1, 'def': 4.},
               'unit': '[R_sun]'},
    'alpha': {'Standard': {'min': 0., 'max': 90., 'step': 0.1, 'def': 45.},
              '<10Rsun': {'min': 0., 'max': 90., 'step': 0.05, 'def': 45.},
              '>10Rsun': {'min': 0., 'max': 90., 'step': 0.25, 'def': 45.},
              '>30Rsun': {'min': 0., 'max': 90., 'step': 0.5, 'def': 45.},
              'unit': '[deg.]'},
    'kappa': {'Standard': {'min': 0., 'max': 1., 'step': 0.01, 'def': 0.3},
              '<10Rsun': {'min': 0., 'max': 1., 'step': 0.01, 'def': 0.3},
              '>10Rsun': {'min': 0., 'max': 1., 'step': 0.01, 'def': 0.3},
              '>30Rsun': {'min': 0., 'max': 1., 'step': 0.01, 'def': 0.3},
              'unit': ''},
    'tilt': {'Standard': {'min': -90., 'max': 90., 'step': 0.25, 'def': 0.},
             '<10Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             '>10Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             '>30Rsun': {'min': -90., 'max': 90., 'step': 0.1, 'def': 0.},
             'unit': '[deg.]'},
}

sliders_dict = {'Spheroid': sliders_Spheroid,
                'Ellipsoid': {**sliders_Spheroid, **sliders_Ellipsoid},
                'GCS': sliders_GCS}


############################################################
# Plus and minus to the images c-limits

slider_image_pmclims = {
    'Running Diff.': [-50, 50],
    'Base Diff.': [-100, 100],
    'Plain': [-5, 5],
}
