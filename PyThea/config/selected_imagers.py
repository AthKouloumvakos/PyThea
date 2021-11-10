'''
 A dictionary of the selected imagers
 How to use this:

'''
import astropy.units as u
from sunpy.net import attrs as a

imager_dict = dict.fromkeys(['LC2'], [(a.Instrument.lasco, a.Detector.c2),
                                      {'dimensions': (1024*u.pixel, 1024*u.pixel), 'polar': 'Clear'}, 'LASCO-C2'])
imager_dict.update(dict.fromkeys(['LC3'], [(a.Instrument.lasco, a.Detector.c3),
                                           {'dimensions': (1024*u.pixel, 1024*u.pixel), 'polar': 'Clear'}, 'LASCO-C3']))
imager_dict.update(dict.fromkeys(['AIA'], [(a.Instrument.aia, a.Wavelength(19.3 * u.nm), a.Sample(2.5*u.minute)),
                                           {'dimensions': (4096*u.pixel, 4096*u.pixel)}, 'SDO-AIA']))
imager_dict.update(dict.fromkeys(['COR2A'], [(a.Source('STEREO_A'), a.Detector.cor2),
                                             {'dimensions': (2048*u.pixel, 2048*u.pixel), 'polar': 1001}, 'STA-COR2']))
imager_dict.update(dict.fromkeys(['COR2B'], [(a.Source('STEREO_B'), a.Detector.cor2),
                                             {'dimensions': (2048*u.pixel, 2048*u.pixel), 'polar': 1001}, 'STB-COR2']))
imager_dict.update(dict.fromkeys(['EUVIA'], [(a.Source('STEREO_A'), a.Detector.euvi, a.Wavelength(19.5 * u.nm)),
                                             {'dimensions': (2048*u.pixel, 2048*u.pixel)}, 'STA-EUVI']))
imager_dict.update(dict.fromkeys(['EUVIB'], [(a.Source('STEREO_B'), a.Detector.euvi, a.Wavelength(19.5 * u.nm)),
                                             {'dimensions': (2048*u.pixel, 2048*u.pixel)}, 'STB-EUVI']))
imager_dict.update(dict.fromkeys(['COR1A'], [(a.Source('STEREO_A'), a.Detector.cor1),
                                             {'dimensions': (512*u.pixel, 512*u.pixel), 'polar': 0}, 'STA-COR1']))
imager_dict.update(dict.fromkeys(['COR1B'], [(a.Source('STEREO_B'), a.Detector.cor1),
                                             {'dimensions': (512*u.pixel, 512*u.pixel), 'polar': 0}, 'STB-COR1']))
imager_dict.update(dict.fromkeys(['HI1A'], [(a.Source('STEREO_A'), a.Detector.hi1),
                                            {'dimensions': (1024*u.pixel, 1024*u.pixel)}, 'STA-HI1']))
imager_dict.update(dict.fromkeys(['HI1B'], [(a.Source('STEREO_B'), a.Detector.hi1),
                                            {'dimensions': (1024*u.pixel, 1024*u.pixel)}, 'STB-HI1']))
imager_dict.update(dict.fromkeys(['HI2A'], [(a.Source('STEREO_A'), a.Detector.hi2),
                                            {'dimensions': (1024*u.pixel, 1024*u.pixel)}, 'STA-HI2']))
imager_dict.update(dict.fromkeys(['HI2B'], [(a.Source('STEREO_B'), a.Detector.hi2),
                                            {'dimensions': (1024*u.pixel, 1024*u.pixel)}, 'STB-HI2']))
