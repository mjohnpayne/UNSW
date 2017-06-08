#! /usr/bin/env python
import os
"""

Downloads all MLST Alleles sequences, and ST profiles from EnteroBase

Includes University of Warwick (formally UCC) MLST Schemes curated by
Mark Achtman.(http://mlst.warwick.ac.uk)

For rMLST please contact University of Oxford directly (http://pubmlst.org/rmlst/). 


### CHANGE LOG ###
2016-01-19 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Initial build
2016-11-17 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Modified to work with api v2.0
    * Fast download of daily dump of database 
    * Incremental download of additional STs and alleles added since.
    * Added better opts to configure 
"""

import os
import logging
import traceback
import time
import argparse
import urllib2
import json
import base64
import sys
from urllib2 import HTTPError
import textwrap
import gzip
from random import shuffle

__epi__ = "Licence: GPLv3 by Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>"
__version__ = '0.1.1'

# You must have a valid API Token
API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)
SERVER_ADDRESS = 'http://enterobase.warwick.ac.uk'


def __create_request(request_str):

    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (API_TOKEN,'')).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    return request

def main():
    global args
    DATABASE = args.database 
    scheme = args.scheme
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if args.listschemes:
        list_schemes(DATABASE)
    else:
        if args.dump:
            dump_scheme(DATABASE, scheme) 
        get_st_profiles(DATABASE, scheme, args.reldate)
        get_alleles(DATABASE, scheme, args.reldate)  
        
        
def dump_scheme(DATABASE, scheme) :
    address = SERVER_ADDRESS + '/api/v2.0/%s/schemes?scheme_name=%s&limit=%d&only_fields=download_st_link' %(DATABASE, scheme, 4000)
    logger.info('Retrieving scheme information from %s' %address)
    try:
        response = urllib2.urlopen(__create_request(address))
        data = json.load(response)
        for scheme_record in data['Schemes']:
            profile_link = scheme_record.get('download_sts_link', None)
            if profile_link:
                logger.info('Downloading %s ST Profiles' %scheme)
                response = urllib2.urlopen(profile_link)
                with open(os.path.join(args.output, '%s-profiles.gz' %scheme), 'wb') as output_profile:
                    output_profile.write(response.read())
        address = SERVER_ADDRESS + '/api/v2.0/%s/%s/loci?limit=%d' %(DATABASE, scheme, 4000)
        logger.info('Retrieving scheme information from %s' %address)
        response = urllib2.urlopen(__create_request(address))
        data = json.load(response)
        for idx, locus in enumerate(data['loci']):
            allele_link = locus.get('download_alleles_link', None)
            if allele_link:
                if idx % 500 == 0 or idx == len(data['loci'])-1:
                    logger.info('Downloading %d out of %d alleles' %(idx, len(data['loci'])))
                response = urllib2.urlopen(allele_link)
                with open(os.path.join(args.output, '%s-%s.gz' %(scheme,str(locus['locus']))), 'wb') as output_allele:
                    output_allele.write(response.read())
                output_allele.close()
    except HTTPError as Response_error:
        logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                  Response_error.msg,
                                                  Response_error.geturl(),
                                                  Response_error.read()))      
        
def list_schemes(DATABASE):
    address = SERVER_ADDRESS + '/api/v2.0/%s/schemes?limit=%d' %(DATABASE, 4000)
    logger.info('Retrieving scheme information from %s' %address)        
    try:
        response = urllib2.urlopen(__create_request(address))
        data = json.load(response)     
        print json.dumps(data, indent=4, sort_keys=True)
        for idx, scheme in enumerate(data['Schemes']):
            if scheme['download_sts_link']:
                import pprint
                print pprint.pprint(scheme)
               # logger.info('Scheme %d: %s\t(%s)' %(idx, scheme['scheme_name'], scheme['label'] ))
    except HTTPError as Response_error:
        logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                  Response_error.msg,
                                                  Response_error.geturl(),
                                                  Response_error.read()))      
def get_alleles(DATABASE, scheme, reldate=1):
    if API_TOKEN:
        try:
            if not os.path.exists(args.output):
                os.mkdir(args.output)
            limit = 5000
            address = SERVER_ADDRESS + '/api/v2.0/%s/%s/loci?only_fields=locus&limit=%d' %(DATABASE, scheme, limit)
            response = urllib2.urlopen(__create_request(address))
            data = json.load(response)
            loci = [] 
            for locus in data['loci']:
                if locus.get('locus'):
                    loci.append(locus.get('locus'))
            chunks = [loci[x:x+10] for x in xrange(0, len(loci), 10)]
            count = 0 
            for locus in chunks:
                offset = 0
                limit = 500
                while True:
                    while True:
                        address = SERVER_ADDRESS + '/api/v2.0/%s/%s/alleles?locus=%s&reldate=%d&only_fields=allele_id,locus,seq&limit=%d&offset=%d' %(DATABASE, scheme, ','.join(locus), reldate, limit,offset)
                        try:
                            response = urllib2.urlopen(__create_request(address))
                            if response.code == 200:
                                break
                        except Exception: 
                            logger.error('Error reading %s' %address)
                    offset += limit
                    data = json.load(response)
                    logger.info('Fetching  %d' %offset)
                    if data['alleles']:
                        _write_allele_file(scheme, data['alleles'])
                    if data['links']['paging'].get('next') == None:
                        break
                count += len(locus)
                logger.info('Loaded %d loci' %count)
            logger.info('Output in directory %s' %args.output)

        except HTTPError as Response_error:
            logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                      Response_error.msg,
                                                      Response_error.geturl(),
                                                      Response_error.read()))
    else:
        logger.error('No API Token, Please register at enterobase.warwick.ac.uk '
                     'and request a token from an administrator')

def _write_allele_file(scheme, alleles):
    groups = {}
    for allele in alleles:
        if groups.get(allele['locus'], None):
            groups[allele['locus']].append(allele)
        else:
            groups[allele['locus']] = [allele]
    for allele_list in groups:
        allele_file_loc = os.path.join(args.output, '%s-%s.gz' %(scheme,allele_list))
        existing_alleles = [] 
        if os.path.exists(allele_file_loc):
            with gzip.open(allele_file_loc, 'r') as allele_reader:
                for line in allele_reader.readlines():
                    if line.startswith('>'):
                        existing_alleles.append(line[1:].strip())
            allele_file = gzip.open(allele_file_loc, 'a')
        else:
            allele_file = gzip.open(allele_file_loc, 'w')                    
        for allele in groups[allele_list]:
            header = '%s_%s' %(allele['locus'], allele['allele_id'])
            if allele['allele_id'] > 0 and header not in existing_alleles:
                    allele_file.write('>%s\n%s\n' %(header,
                                             textwrap.fill(allele['seq'], 80)))
        allele_file.close()

def get_st_profiles(DATABASE, scheme, reldate=1):
    if API_TOKEN:        
        try:
            st_file_loc = os.path.join(args.output, '%s-profiles.gz' %scheme)
            existing_sts = []
            if os.path.exists(st_file_loc):
                with gzip.open(st_file_loc, 'r') as st_reader:
                    for line in st_reader.readlines()[1:]:
                        existing_sts.append(line.split('\t')[0])
                st_file = gzip.open(st_file_loc, 'a')
                first = False
            else:
                st_file = gzip.open(st_file_loc, 'w')
                first = True
            offset = 0
            limit = 700
            while True:
                address = SERVER_ADDRESS + '/api/v2.0/%s/%s/sts?show_alleles=true&reldate=%d&limit=%d&offset=%d' %(DATABASE, scheme, reldate, limit, offset)
                response = urllib2.urlopen(__create_request(address))
                data = json.load(response)
                offset += limit
                for st in data['STs']:
                    if st['ST_id'] > 0 :
                        allele_dict = {}
                        for allele in st['alleles']:
                            allele_dict[allele['locus']] = allele['allele_id']
                        allele_row = []
                        for key in sorted(allele_dict.keys()):
                            allele_row.append(str(allele_dict[key]))
                        if first:
                            st_file.write('ST_id\t%s\n' %'\t'.join(sorted(allele_dict.keys())))
                            first = False
                        if st['ST_id'] not in existing_sts:
                            st_file.write('%s\t%s\n' %(st['ST_id'],
                                                           '\t'.join(allele_row)))

                st_file.flush()
                logger.info('Loaded %d of %d'%(offset, data['links']['total_records']))
                if data['links']['paging'].get('next') == None:
                    break
            st_file.close()
            logger.info('Output in directory %s' %args.output)
        except HTTPError as Response_error:
            logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                      Response_error.msg,
                                                      Response_error.geturl(),
                                                      Response_error.read()))
    else:
        logger.error('No API Token, Please register at enterobase.warwick.ac.uk '
                     'and request a token from an administrator')
    


if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = "cgMLST_mp" #__doc__.split('\n\n')[1].strip()
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        parser = argparse.ArgumentParser(description=desc,epilog=__epi__)
        parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        parser.add_argument('-d','--database',action='store', help='Database name either senterica,ecoli,yersinia', 
                            choices=['senterica','ecoli','yersinia'], default='senterica')
        parser.add_argument('-s','--scheme',action='store', help='Scheme name to download e.g. cgMLST, MLST_Achtman', default='MLST_Achman')
        parser.add_argument('-x','--dump',action='store_true', help='Download daily dump of all records', default=False)        
        parser.add_argument('-r','--reldate',action='store', type=int, help='Retrieve records from the last N days', default=1)        
        parser.add_argument('-l','--listschemes',action='store_true', help='List schemes to choose from and exist')
        parser.add_argument('-o','--output',action='store',help='Output directory', default='entero_out')
        
        args = parser.parse_args()
        if args.verbose: print "Executing @ " + time.asctime()
        main()


        if args.verbose: print "Ended @ " + time.asctime()
        if args.verbose: print 'total time in minutes:',
        if args.verbose: print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)
