from collections import OrderedDict

import pandas as pd
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.time import TimeRange

__all__ = ['STEREOClient']


class STEREOClient(GenericClient):
    """
    Provides access to data from NASA/NASCOM for SECCHI instrument suite on board STEREO.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> res = Fido.search(a.Time('2015-06-21 00:00', '2015-06-23 23:59'),
    ...                    a.Detector('cor2'),
                           a.Instrument('STEREO'),
                           a.Provider('NASCOM'))

    """

    baseurl = r'https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/a/img/INSTRUMENT/%Y%m%d/'

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
        summaryurl = r'https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/??/summary/SCC%Y%m.img.DETECTOR'

        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])

        metalist = []
        print(matchdict['Source'][0])
        if matchdict['Source'][0].upper() == 'STEREO_A':
            scl, scc = 'a', 'sccA'
        elif matchdict['Source'][0].upper() == 'STEREO_B':
            scl, scc = 'b', 'sccB'

        if matchdict['Detector'][0].upper() == 'COR2':
            detect = 'c2'

        for times in pd.date_range(tr.start.strftime('%Y-%m-%d %H:%M:%S'), tr.end.strftime('%Y-%m-%d %H:%M:%S')):
            urlsum =  summaryurl.replace('??', scl).replace('SCC', scc).replace('DETECTOR', detect)
            urlsum = times.strftime(urlsum)  # sccA201105.img.h1

            data = pd.read_fwf(urlsum,
                               header=None, colspecs=([0, 25], [28, 47], [50, 54]), skiprows=2,
                               names=['filename', 'date', 'telescope'])
            data = data.set_index([pd.to_datetime(data.iloc[:, 1], format='%Y/%m/%d %H:%M:%S')])
            data = data[(data.index >= tr.start.strftime('%Y-%m-%d %H:%M:%S')) &
                        (data.index <= tr.end.strftime('%Y-%m-%d %H:%M:%S'))]

            # URL for example: files = 'https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/a/img/cor2/20211028/' + files
            urldata = times.strftime(baseurl).replace('INSTRUMENT', matchdict['Detector'][0])

            for f in data.itertuples():
                filemeta = {'start': f.Index,
                            'url': urldata + f.filename}
                rowdict = self.post_search_hook(filemeta, matchdict)
                metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {
            attrs.Instrument: [('SECCHI', '')],
            attrs.Source: [('STEREO_A', ''), ('STEREO_B', '')],
            attrs.Provider: [('NASCOM', '')],
            attrs.Detector: [('COR2', ''), ('COR1', '')]
        }
        return adict
