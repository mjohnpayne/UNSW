from time import sleep as sl

from urllib2 import HTTPError
import urllib2
import base64
import json
import os

API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)

print API_TOKEN

def __create_request(request_str):

    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (API_TOKEN,'')).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    return request

address = 'http://enterobase.warwick.ac.uk/api/v2.0/senterica/straindata?limit=1&assembly_status=Assembled'

try:
    response = urllib2.urlopen(__create_request(address))
    data = json.load(response)
    print json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))

except HTTPError as Response_error:
    print '%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                      Response_error.msg,
                                                      Response_error.geturl(),
                                                      Response_error.read())