#! /bin/env python3

# packages
import os, sys, re
import numpy as np
import pandas as pd
import sqlalchemy as sa
import urllib.request as url
from Bio import SeqIO
from Bio import Seq
from io import StringIO
# ------------------
from biofuncs import gfftable


# ------------------------------------------------------------------------------
# Class

# ------------------
# UniAPI
class UniAPI:
    # {{{

    '''
UniAPI is used to abtain information from uniprot website using uniprot API.

UniAPI is the basic of several other classes.
'''
    def __init__(self, part = 'uniprot'):
        self._web = 'http://www.uniprot.org/'
        if part not in ['uniprot', 'uniparc', 'uniref', 'taxonomy',
                        'locations', 'diseases']:
            raise ValueError(
                'Wrong value for "part" parameter: {0:s}. \n("part" should be in ["uniprot", "uniparc", "uniref"])'.format(part)
            )
        self._web = self._web + part

        # }}}


# ------------------
# UniProt
class UniProt(UniAPI):
    # {{{

    '''
UniProt is used to abtain protein information by uniprot id of proteins.'''
    def __init__(self, uniprot_id):
        UniAPI.__init__(self, part = 'uniprot')
        self.uniprot_id = uniprot_id
        self.__sequence()
        self.__anno()
        self.detail()
    def _get_info(self, format):
        '''Fetch data of protein from uniprot by uniprot id. Result might be xml, txt, fasta, gff.'''
        if format not in ['xml', 'txt', 'fasta', 'gff']:
            raise ValueError(
                'Wrong value for "format" parameter: {0:s}. ("format" should be in ["xml", "txt", "fasta", "gff"])'.format(part)
            )
        else:
            proturl = self._web + self._part + '/' + self.uniprot_id + '.' + format
            prot = url.urlopen(urlparse.quote(proturl, safe = '/:+=,&?'))
            if prot.getcode() == 200:
                if format == 'fasta':
                    self.seq = SeqIO.read(
                        StringIO(prot.read().decode('utf-8')),
                        'fasta'
                    )
                elif format == 'gff':
                    gffcolumn = [
                        'uniprot_id',     'source', 'feature',
                        'start',  'end',    'score',
                        'strand', 'frame',  'attribute'
                    ]
                    self.gff = gff_table(
                        proturl,
                        'gff',
                        gffcolumn
                    )
                elif format == 'xml':
                    self.xml = et.parse(StringIO(prot.read().decode('utf-8')))                    
                elif format == 'txt':
                    self.txt = prot.read().decode('utf-8')
                else:
                    pass
            elif prot.getcode() == 404:
                raise ValueError(
                    'Not found. The resource you requested doesn’t exist.'
                )
            elif prot.getcode() == 410:
                raise ValueError(
                    'Gone. The resource you requested was removed.'
                )
            elif prot.getcode() == 500:
                raise ValueError(
                    'Internal server error. Most likely a temporary problem, but if the problem persists please contact uniprot.'
                )
            elif prot.getcode() == 503:
                raise ValueError(
                    'Service not available. The server is being updated, try again later.'
                )
            else:
                raise ValueError(
                    'Bad request. There is a problem with your input.'
                )
    def __sequence(self):
        if not hasattr(self, 'seq'):
            self._get_info('fasta')
        return self.seq
    def __anno(self):
        if not hasattr(self, 'gff'):
            self._get_info('gff')
        return self.gff
    def detail(self, raw = False):
        if not hasattr(self, 'xml'):
            self._get_info('xml')
        elif not hasattr(self, 'txt'):
            self._get_info('txt')
        if raw:
            return self.txt
        else:
            return self.xml
    def __getattr__(self, attr):
        if attr == 'sequence':
            return self.seq
        elif attr == 'annotation' or attr == 'anno':
            return self.gff
        elif attr == 'infomation':
            return self.xml
        else:
            raise AttributeError(attr)

    # }}}


# ------------------
# UniQuery
class UniQuery(UniAPI):
    # {{{

    '''
UniQuery is used to query from uniprot by url. The query string must be appropriate.
'''
    uniprot_format = [
        'html', 'tab', 'xls', 'fasta', 'gff',
        'txt', 'xml', 'rdf', 'list', 'rss'
    ]
    uniprot_column = [
        'id',
        'entry name',
        'genes',
        'genes(PREFERRED)',
        'genes(ALTERNATIVE)',
        'genes(OLN)',
        'genes(ORF)',
        'organism',
        'organism-id',
        'protein names',
        'proteome',
        'lineage(ALL)',
        'virus hosts',
        'fragment',
        'encodedon',
        'comment(ALTERNATIVE PRODUCTS)',
        'comment(ERRONEOUS GENE MODEL PREDICTION)',
        'comment(ERRONEOUS INITIATION)',
        'comment(ERRONEOUS TERMINATION)',
        'comment(ERRONEOUS TRANSLATION)',
        'comment(FRAMESHIFT)',
        'comment(MASS SPECTROMETRY)',
        'comment(POLYMORPHISM)',
        'comment(RNA EDITING)',
        'comment(SEQUENCE CAUTION)',
        'length',
        'mass',
        'sequence',
        'feature(ALTERNATIVE SEQUENCE)',
        'feature(NATURAL VARIANT)',
        'feature(NON ADJACENT RESIDUES)',
        'feature(NON STANDARD RESIDUE)',
        'feature(NON TERMINAL RESIDUE)',
        'feature(SEQUENCE CONFLICT)',
        'feature(SEQUENCE UNCERTAINTY)',
        'version(sequence)',
        'ec',
        'comment(ABSORPTION)',
        'comment(CATALYTIC ACTIVITY)',
        'comment(COFACTOR)',
        'comment(ENZYME REGULATION)',
        'comment(FUNCTION)',
        'comment(KINETICS)',
        'comment(PATHWAY)',
        'comment(REDOX POTENTIAL)',
        'comment(TEMPERATURE DEPENDENCE)',
        'comment(PH DEPENDENCE)',
        'feature(ACTIVE SITE)',
        'feature(BINDING SITE)',
        'feature(DNA BINDING)',
        'feature(METAL BINDING)',
        'feature(NP BIND)',
        'feature(SITE)',
        'annotation score',
        'features',
        'comment(CAUTION)',
        'comment(GENERAL)',
        'keywords',
        'context',
        'existence',
        'tools',
        'reviewed',
        'comment(SUBUNIT)',
        'interactor',
        'comment(DEVELOPMENTAL STAGE)',
        'comment(INDUCTION)',
        'comment(TISSUE SPECIFICITY)',
        'go',
        'go(biological process)',
        'go(molecular function)',
        'go(cellular component)',
        'go-id',
        'comment(ALLERGEN)',
        'comment(BIOTECHNOLOGY)',
        'comment(DISRUPTION PHENOTYPE)',
        'comment(DISEASE)',
        'comment(PHARMACEUTICAL)',
        'comment(TOXIC DOSE)',
        'comment(PTM)',
        'feature(CHAIN)',
        'feature(CROSS LINK)',
        'feature(DISULFIDE BOND)',
        'feature(GLYCOSYLATION)',
        'feature(INITIATOR METHIONINE)',
        'feature(LIPIDATION)',
        'feature(MODIFIED RESIDUE)',
        'feature(PEPTIDE)',
        'feature(PROPEPTIDE)',
        'feature(SIGNAL)',
        'feature(TRANSIT)',
        '3d',
        'feature(BETA STRAND)',
        'feature(HELIX)',
        'feature(TURN)',
        'citationmapping',
        'citation',
        'created',
        'last-modified',
        'sequence-modified',
        'version(entry)',
        'comment(DOMAIN)',
        'comment(SIMILARITY)',
        'families',
        'feature(COILED COIL)',
        'feature(COMPOSITIONAL BIAS)',
        'feature(DOMAIN EXTENT)',
        'feature(MOTIF)',
        'feature(REGION)',
        'feature(REPEAT)',
        'feature(ZINC FINGER)',
        'lineage(all)',
        'lineage(SUPERKINGDOM)',
        'lineage(KINGDOM)',
        'lineage(SUBKINGDOM)',
        'lineage(SUPERPHYLUM)',
        'lineage(PHYLUM)',
        'lineage(SUBPHYLUM)',
        'lineage(SUPERCLASS)',
        'lineage(CLASS)',
        'lineage(SUBCLASS)',
        'lineage(INFRACLASS)',
        'lineage(SUPERORDER)',
        'lineage(ORDER)',
        'lineage(SUBORDER)',
        'lineage(INFRAORDER)',
        'lineage(PARVORDER)',
        'lineage(SUPERFAMILY)',
        'lineage(FAMILY)',
        'lineage(SUBFAMILY)',
        'lineage(TRIBE)',
        'lineage(SUBTRIBE)',
        'lineage(GENUS)',
        'lineage(SUBGENUS)',
        'lineage(SPECIES GROUP)',
        'lineage(SPECIES SUBGROUP)',
        'lineage(SPECIES)',
        'lineage(SUBSPECIES)',
        'lineage(VARIETAS)',
        'lineage(FORMA)',
        'lineage-id(all)',
        'lineage-id(SUPERKINGDOM)',
        'lineage-id(KINGDOM)',
        'lineage-id(SUBKINGDOM)',
        'lineage-id(SUPERPHYLUM)',
        'lineage-id(PHYLUM)',
        'lineage-id(SUBPHYLUM)',
        'lineage-id(SUPERCLASS)',
        'lineage-id(CLASS)',
        'lineage-id(SUBCLASS)',
        'lineage-id(INFRACLASS)',
        'lineage-id(SUPERORDER)',
        'lineage-id(ORDER)',
        'lineage-id(SUBORDER)',
        'lineage-id(INFRAORDER)',
        'lineage-id(PARVORDER)',
        'lineage-id(SUPERFAMILY)',
        'lineage-id(FAMILY)',
        'lineage-id(SUBFAMILY)',
        'lineage-id(TRIBE)',
        'lineage-id(SUBTRIBE)',
        'lineage-id(GENUS)',
        'lineage-id(SUBGENUS)',
        'lineage-id(SPECIES GROUP)',
        'lineage-id(SPECIES SUBGROUP)',
        'lineage-id(SPECIES)',
        'lineage-id(SUBSPECIES)',
        'lineage-id(VARIETAS)',
        'lineage-id(FORMA)'
    ]
    def __init__(self,
                 part     = 'uniprot',
                 query    = '',
                 format   = 'tab',
                 columns  = ['id', 'entry name', 'protein names', 'genes(PREFERRED)', 'organism'],
                 include  = 'yes',
                 compress = 'no',
                 limit    = None,
                 offset   = None):
        self._query    = query
        self._format   = format
        self._columns  = columns
        self._include  = include
        self._compress = compress
        self._limit    = limit
        self._offset   = offset
        UniAPI.__init__(self, part)
        if self._format not in UniQuery.uniprot_format: # format have invalid values.
            raise ValueError(
                'Wrong value for "format" parameter: {0:s}. \n ("format" should be in ["html", "tab", "xls", "fasta", "gff", "txt", "xml", "rdf", "list", "rss"])'.format(self._format)
            )
        if sum(x in UniQuery.uniprot_column for x in self._columns) < len(self._columns): # if columns have invalid values.
            raise ValueError(
                'Wrong value for "columns" parameter. \n(Please check http://www.uniprot.org/help/uniprotkb_column_names)'
            )
        if self._include not in ['yes', 'no']:
            raise ValueError(
                'Wrong value for "include" parameter: {0:s}. \n(choose from "yes" or "no")'.format(self._include)
            )
        if self._compress not in ['yes', 'no']:
            raise ValueError(
                'Wrong value for "compress" parameter: {0:s}. \n(choose from "yes" or "no")'.format(self._compress)
            )
        self.__make_url()
        self.__query()

    def __repr__(self):
        return '<uniprot.UniQuery: {0:s}>'.format(self._query)

    def __str__(self):
        return '\n'.join(
            ['UniQuery',
             '--------------------',
             'part: {0:s}',
             'query: {1:s}',
             'format: {2:s}',
             'columns: {3:s}',
             'url: {4:s}']
        ).format(
            self._web,
            self._query,
            self._format,
            ','.join(list(self._columns)),
            urlparse.quote(self._url, safe = '/:+=,&?')
        )
        
    def __make_url(self):
        '''
Make the whole url to query. Without change into percent encode.
'''
        self._url = ''.join(
            [self._web,
             '/?query=', self._query if self._query != '' else '*',
             '&format=', self._format,
             '&include=', self._include,
             '&compress=', self._compress,
             '&limit=', str(self._limit),
             '&offset=', str(self._offset),
             '&sort=score']
        )
        if self._columns != '':
            self._url = '&columns='.join(
                [self._url, ','.join(list(self._columns))]
            )
        return self._url
    
    def reset_query(self, query, **kargs):
        '''
Modify query.
'''
        self._query = '+'.join(
            [query,
             '+'.join(
                 ':'.join([key, value.join(['(', ')'])]) for key, value in kargs.items()
             )]
        )
        self.__make_url()
        self.__query()
        return self._url
    
    def __query(self):
        '''
Query UniProt web.
'''
        prot = url.urlopen(urlparse.quote(self._url, safe = '/:+=,&?'))
        if not hasattr(self, 'result'):
            if prot.getcode() == 200:
                self.result = pd.read_table(prot)
        return self.result

    # }}}


# ------------------
# UniSets
class UniSets():
    # {{{

    def __init__(self):
        self._uni_url = {
            'location': UniQuery(
                part     = 'locations',
                query    = '',
                format   = 'tab',
                columns  = ['id'],
                include  = 'yes',
                compress = 'no',
                limit    = None,
                offset   = None
            ),
            'disease': UniQuery(
                part     = 'diseases',
                query    = '',
                format   = 'tab',
                columns  = ['id'],
                include  = 'yes',
                compress = 'no',
                limit    = None,
                offset   = None
            ),
            'taxonomy': UniQuery(
                part     = 'taxonomy',
                query    = 'reviewed:yes',
                format   = 'tab',
                columns  = ['id'],
                include  = 'yes',
                compress = 'no',
                limit    = None,
                offset   = None
            )
        }
    def list(self):
        '''List posible sets.'''
        return list(self._uni_url.keys())
    def __getitem__(self, setname):
        '''Return the DataFrame of set by setname.'''
        if setname not in self._uni_url:
            raise ValueError(
                'Wrong value for "setname" parameter: {0:s}. \n(choose from [{1:s}])'.format(
                    setname,
                    ', '.join(x for x in self._uni_url.keys())
                )
            )            
        return self._uni_url[setname]
    def __call__(self, setname):
        return self.__getitem__(setname).result

    # }}}


# ------------------
# UniLocation
class UniLocation(UniQuery):
    # {{{

    def __init__(self,
                 location_id   = '',
                 location_name = '',
                 organism      = '',
                 reviewed      = 'yes',
                 columns       = ['id', 'entry name', 'protein names',
                                  'genes(PREFERRED)', 'organism',
                                  'comment(SUBCELLULAR LOCATION)'],
                 limit         = 10000,
                 offset        = 0):
        if location_id != '' and location_name != '':
            self._location = ''.join(
                ['location:"', location_name.capitalize(), ' [',
                 location_id.upper(), ']"']
            )
        elif location_id != '' or location_name != '':
            self._location = ''.join(
                ['location:',
                 location_name.capitalize() if location_name != '' else location_id.upper()]
            )
        else:
            raise ValueError(
                'Please enter location_id or location_name.'
            )
        self._organism = organism
        UniQuery.__init__(
            part     = 'uniprot',
            query    = '',
            format   = 'tab',
            columns  = columns,              
            include  = 'yes',
            compress = 'no',
            limit    = limit,
            offset   = offset
        )
        reviewed = 'yes' if reviewed == 'yes' else 'no'
        self.reset_query(
            ''.join(['locations:(', self._location, ')' ]),
             reviewed = reviewed,
             organism = organism
        )

    def __repr__(self):
        return '\n'.join(
            ['<UniLocation: {0:s}>'.format(self._location),
             '{0:d} records in total.'.format(self.result.shape[0])]
        )

    def __str__(self):
        return self.__repr__()

    def __call__(self):
        return self.result

    # }}}
        
# ------------------------------------------------------------------------------







