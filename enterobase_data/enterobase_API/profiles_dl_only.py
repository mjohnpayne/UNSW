from time import sleep as sl
from urllib2 import HTTPError
import urllib2
import base64
import json
import os

SERVER_ADDRESS = 'http://enterobase.warwick.ac.uk'
DATABASE = 'senterica'
scheme = 'cgMLSTv2'
API_TOKEN = "eyJhbGciOiJIUzI1NiIsImV4cCI6MTUxMjY2MDMxOSwiaWF0IjoxNDk2ODkyMzE5fQ.eyJ1c2VybmFtZSI6Ik1qb2hucGF5bmUiLCJjaXR5IjoiU3lkbmV5IiwiY29uZmlybWVkIjoxLCJhbGxvd2VkX3NjaGVtZXMiOiJjZ01MU1RfdjIiLCJmaXJzdG5hbWUiOiJNaWNoYWVsIiwiYXBpX2FjY2Vzc19zZW50ZXJpY2EiOiJUcnVlIiwiY291bnRyeSI6IkF1c3RyYWxpYSIsImlkIjoxMTA1LCJhZG1pbmlzdHJhdG9yIjpudWxsLCJlbWFpbCI6Im1pY2hhZWwucGF5bmVAdW5zdy5lZHUuYXUiLCJkZXBhcnRtZW50IjoiQmlvdGVjaG5vbG9neSBhbmQgYmlvbW9sZWN1bGFyIHNjaWVuY2VzIiwidmlld19zcGVjaWVzIjoiVHJ1ZSIsImxhc3RuYW1lIjoiUGF5bmUiLCJhY3RpdmUiOm51bGwsInVwbG9hZF9yZWFkcyI6IlRydWUiLCJpbnN0aXR1dGlvbiI6IlRoZSBVbml2ZXJzaXR5IG9mIE5ldyBTb3V0aCBXYWxlcyJ9.WLDMlYWJn61zx7dKY-JVAzeM-3CPZyAaYqM74B3kNLo"

def __create_request(request_str):

    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (API_TOKEN,'')).replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    return request

address = SERVER_ADDRESS + '/api/v2.0/%s/schemes?scheme_name=%s&limit=%d&only_fields=download_sts_link' %(DATABASE, scheme, 4000)

os.mkdir(scheme)
try:
    response = urllib2.urlopen(__create_request(address))
    data = json.load(response)
    for scheme_record in data['Schemes']:
        profile_link = scheme_record.get('download_sts_link', None)
        if profile_link:
           response = urllib2.urlopen(profile_link)
           with open(os.path.join(scheme, 'snps-profiles.gz'), 'wb') as output_profile:
               output_profile.write(response.read())
except HTTPError as Response_error:
    print '%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                      Response_error.msg,
                                                      Response_error.geturl(),
                                                      Response_error.read())