from collections import OrderedDict

import pandas as pd
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.time import TimeRange

__all__ = ['LASCOClient']


class LASCOClient(GenericClient):
    """
    Provides access to data from NASA/NASCOM for LASCO instrument on board SOHO.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> res = Fido.search(a.Time('2015-06-21 00:00', '2015-06-23 23:59'),
    ...                    a.Detector('c2'),
                           a.Instrument('LASCO'),
                           a.Provider('NASCOM'))

    """

    baseurl = r'https://sohoftp.nascom.nasa.gov/qkl/lasco/quicklook/level_05/%y%m%d/'
    # r'https://umbra.nascom.nasa.gov/pub/lasco_level05/%y%m%d/'

    def post_search_hook(self, i, matchdict):

        rowdict = OrderedDict()
        rowdict['Start Time'] = i['start']
        rowdict['End Time'] = i['start']
        rowdict['Instrument'] = matchdict['Instrument'][0].upper()
        rowdict['Source'] = matchdict['Source'][0]
        rowdict['Provider'] = matchdict['Provider'][0]
        rowdict['url'] = i['url']

        return rowdict

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
             Any extra keywords to refine the search.

        Returns
        -------
        A `QueryResponse` instance containing the query result.
        """

        baseurl, _, matchdict = self.pre_search_hook(*args, **kwargs)

        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])

        metalist = []

        for times in pd.date_range(tr.start.strftime('%Y-%m-%d %H:%M:%S'), tr.end.strftime('%Y-%m-%d %H:%M:%S')):
            url = times.strftime(baseurl)

            data = pd.read_fwf(url + matchdict['Detector'][0] + '/img_hdr.txt',
                               header=None, widths=(13, 13, 13),
                               names=['filename', 'date', 'time'])
            data = data.set_index([pd.to_datetime(data.iloc[:, 1]+' '+data.iloc[:, 2], format='%Y/%m/%d %H:%M:%S')])
            data = data[(data.index >= tr.start.strftime('%Y-%m-%d %H:%M:%S')) &
                        (data.index <= tr.end.strftime('%Y-%m-%d %H:%M:%S'))]
            # files = 'https://sohoftp.nascom.nasa.gov/qkl/lasco/quicklook/level_05/210804/' + 'c2/' + files

            for f in data.itertuples():
                filemeta = {'start': f.Index,
                            'url': url+matchdict['Detector'][0]+'/' + f.filename}
                rowdict = self.post_search_hook(filemeta, matchdict)
                metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('LASCO', '')],
                 attrs.Source: [('SOHO', '')],
                 attrs.Provider: [('NASCOM', '')],
                 attrs.Detector: [('c2', ''), ('c3', '')]}
        return adict
