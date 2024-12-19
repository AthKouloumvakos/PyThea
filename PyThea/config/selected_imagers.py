'''
 A dictionary of the selected imagers
 How to use this:

'''
import astropy.units as u
import sunpy_soar  # noqa
from sunpy.net import attrs as a

imager_dict = {}

imager_dict['LC2'] = {'fido': (a.Instrument.lasco, a.Detector.c2),
                      'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'polar': 'Clear', 'superpixel': 2},
                      'source': 'SOHO', 'instrument': 'LASCO', 'detector': 'C2'}

imager_dict['LC3'] = {'fido': (a.Instrument.lasco, a.Detector.c3),
                      'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'polar': 'Clear', 'superpixel': 2},
                      'source': 'SOHO', 'instrument': 'LASCO', 'detector': 'C3'}

imager_dict['AIA-193'] = {'fido': (a.Instrument.aia, a.Wavelength(19.3 * u.nm), a.Sample(1*u.minute)),
                          'process': {'dimensions': (4096*u.pixel, 4096*u.pixel), 'superpixel': 8, 'exposure': 1.90},
                          'source': 'SDO', 'instrument': 'AIA', 'wavelength': '193'}

imager_dict['AIA-211'] = {'fido': (a.Instrument.aia, a.Wavelength(21.1 * u.nm), a.Sample(1*u.minute)),
                          'process': {'dimensions': (4096*u.pixel, 4096*u.pixel), 'superpixel': 8, 'exposure': 1.90},
                          'source': 'SDO', 'instrument': 'AIA', 'wavelength': '211'}

imager_dict['COR2A'] = {'fido': (a.Source('STEREO_A'), a.Detector.cor2),
                        'process': {'dimensions': (2048*u.pixel, 2048*u.pixel), 'polar': 1001, 'superpixel': 4},
                        'source': 'STEREO_A', 'instrument': 'SECCHI', 'detector': 'COR2'}

imager_dict['COR2B'] = {'fido': (a.Source('STEREO_B'), a.Detector.cor2),
                        'process': {'dimensions': (2048*u.pixel, 2048*u.pixel), 'polar': 1001, 'superpixel': 4},
                        'source': 'STEREO_B', 'instrument': 'SECCHI', 'detector': 'COR2'}

imager_dict['EUVIA'] = {'fido': (a.Source('STEREO_A'), a.Detector.euvi, a.Wavelength(19.5 * u.nm)),
                        'process': {'dimensions': (2048*u.pixel, 2048*u.pixel), 'superpixel': 4},
                        'source': 'STEREO_A', 'instrument': 'SECCHI', 'detector': 'EUVI'}

imager_dict['EUVIB'] = {'fido': (a.Source('STEREO_B'), a.Detector.euvi, a.Wavelength(19.5 * u.nm)),
                        'process': {'dimensions': (2048*u.pixel, 2048*u.pixel), 'superpixel': 4},
                        'source': 'STEREO_B', 'instrument': 'SECCHI', 'detector': 'EUVI'}

imager_dict['COR1A'] = {'fido': (a.Source('STEREO_A'), a.Detector.cor1),
                        'process': {'dimensions': (512*u.pixel, 512*u.pixel)},
                        'source': 'STEREO_A', 'instrument': 'SECCHI', 'detector': 'COR1'}

imager_dict['COR1B'] = {'fido': (a.Source('STEREO_B'), a.Detector.cor1),
                        'process': {'dimensions': (512*u.pixel, 512*u.pixel)},
                        'source': 'STEREO_B', 'instrument': 'SECCHI', 'detector': 'COR1'}

imager_dict['HI1A'] = {'fido': (a.Source('STEREO_A'), a.Detector.hi1),
                       'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'superpixel': 2},
                       'source': 'STEREO_A', 'instrument': 'SECCHI', 'detector': 'HI1'}

imager_dict['HI1B'] = {'fido': (a.Source('STEREO_B'), a.Detector.hi1),
                       'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'superpixel': 2},
                       'source': 'STEREO_B', 'instrument': 'SECCHI', 'detector': 'HI1'}

imager_dict['HI2A'] = {'fido': (a.Source('STEREO_A'), a.Detector.hi2),
                       'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'superpixel': 2},
                       'source': 'STEREO_A', 'instrument': 'SECCHI', 'detector': 'HI2'}

imager_dict['HI2B'] = {'fido': (a.Source('STEREO_B'), a.Detector.hi2),
                       'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'superpixel': 2},
                       'source': 'STEREO_B', 'instrument': 'SECCHI', 'detector': 'HI2'}

imager_dict['WISPR1'] = {'fido': (a.Instrument.wispr, a.Detector.inner),
                         'process': {'dimensions': (960*u.pixel, 1024*u.pixel), 'processing_level': 3, 'superpixel': 2},
                         'source': 'PSP', 'instrument': 'WISPR', 'detector': 'Inner'}

imager_dict['WISPR2'] = {'fido': (a.Instrument.wispr, a.Detector.outer),
                         'process': {'dimensions': (960*u.pixel, 1024*u.pixel), 'processing_level': 3, 'superpixel': 2},
                         'source': 'PSP', 'instrument': 'WISPR', 'detector': 'Outer'}

imager_dict['EUI-FSI'] = {'fido': (a.Instrument('EUI'), a.soar.Product('EUI-FSI174-IMAGE'), a.Level(2)),
                          'process': {'superpixel': 1},
                          'source': 'SOLO', 'instrument': 'EUI-FSI', 'wavelength': '174'}

imager_dict['METIS'] = {'fido': (a.Instrument('METIS'), a.soar.Product('METIS-VL-TB'), a.Level(2)),
                        'process': {'dimensions': (1024*u.pixel, 1024*u.pixel), 'superpixel': 2},
                        'source': 'SOLO', 'instrument': 'METIS', 'detector': 'VLD', 'wavelength': 'TB'}

for tile in range(1, 5):
    z = 'T' if tile in [1, 2] else 'G'
    imager_dict[f'SOLOHI-T{tile}'] = {'fido': (a.Instrument('SOLOHI'),
                                               a.soar.Product(f'SOLOHI-{tile}F{z}'), a.Level(2)),
                                      'process': {'superpixel': 2},
                                      'source': 'SOLO', 'instrument': 'SOLOHI', 'detector': f'T{tile}'}
