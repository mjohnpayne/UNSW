#! /usr/bin/env python
"""
Tests your EnteroBase API token with a simple request.

### CHANGE LOG ###
2016-01-19 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Initial build
2016-08-30 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Updated to API 2.0
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

__epi__ = "Licence: GPLv3 by Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>"
__version__ = '0.0.1'

# You must have a valid API Token
API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)
SERVER_ADDRESS = 'http://enterobase.warwick.ac.uk'

def __create_request(request_str):

    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (API_TOKEN,'')).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    return request

def test_api():
    if API_TOKEN:
        global args
        address = SERVER_ADDRESS + '/api/v2.0/senterica/straindata?assembly_status=Assembled&limit=1&offset=10000'
        try:
            response = urllib2.urlopen(__create_request(address))
            data = json.load(response)
            print json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))

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
        desc = __doc__.split('\n\n')[1].strip()
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        parser = argparse.ArgumentParser(description=desc,epilog=__epi__)
        parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        args = parser.parse_args()
        if args.verbose: print "Executing @ " + time.asctime()

        test_api()

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
