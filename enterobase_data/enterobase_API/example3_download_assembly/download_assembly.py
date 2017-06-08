#! /usr/bin/env python
"""

Downloads assemblies from Enterobase

Example 1 uses search parameters to pick strains from a certain serovar. 
Example 2 (-x) downloads all the assemblies in a database

### CHANGE LOG ###
2017-01-16 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Initial build
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

SERVER_ADDRESS = 'http://hydra.warwick.ac.uk:5000'

def __create_request(request_str):

    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (API_TOKEN,'')).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    return request

def main():
    global args
    DATABASE = args.database 
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    address = SERVER_ADDRESS + '/api/v2.0/%s/straindata?serotype=%s&assembly_status=Assembled&limit=%d&only_fields=strain_name,download_fasta_link' %(DATABASE, args.serotype, 40)
    logger.info('Retrieving assembly information from %s' %address)
    try:
        response = urllib2.urlopen(__create_request(address))
        data = json.load(response)
        for record in data['straindata']:
            record_values = data['straindata'][record]
            response = urllib2.urlopen(__create_request(record_values['download_fasta_link']))
            with open(os.path.join(args.output, '%s.fasta' %record_values['strain_name']),'w') as out_ass: 
                out_ass.write(response.read())
    except HTTPError as Response_error:
        logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                  Response_error.msg,
                                                  Response_error.geturl(),
                                                  Response_error.read()))
    if args.dump:
        dump_whole_database(DATABASE)
        
def dump_whole_database(DATABASE):
    address = SERVER_ADDRESS + '/api/v2.0/%s/assemblies?assembly_status=Assembled&limit=%d&only_fields=download_fasta_link,assembly_barcode' %(DATABASE, 1000)
    logger.info('Retrieving assembly information from %s' %address)
    while True:
        try:
            response = urllib2.urlopen(__create_request(address))
            data = json.load(response)
            for record in data['Assemblies']:
                response = urllib2.urlopen(__create_request(record['download_fasta_link']))
                if response.getcode() == 200:
                    with open(os.path.join(args.output, '%s.fasta' %record['assembly_barcode']),'w') as out_ass: 
                        out_ass.write(response.read())
            if data['links']['paging'].get('next'):
                address = data['links']['paging']['next']
            else: 
                break
        except HTTPError as Response_error:
            logger.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                      Response_error.msg,
                                                      Response_error.geturl(),
                                                      Response_error.read()))

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        parser = argparse.ArgumentParser(description=desc,epilog=__epi__)
        parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        parser.add_argument('-d','--database',action='store', help='Database name either senterica,ecoli,yersinia', 
                            choices=['senterica','ecoli'], default='senterica')
        parser.add_argument('-o','--output',action='store',help='Output directory', default='entero_out')
        parser.add_argument('-s','--serotype',action='store',help='Serotype', default='Agona')
        parser.add_argument('-x','--dump',action='store_true',help='Dump whole database', default=False)
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
