#! /usr/bin/env python
"""
Helper Daemon for Bionumerics

### CHANGE LOG ###
2016-08-30 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Initial build

"""

import os
import subprocess
import time
import traceback
import logging
import argparse
import datetime
import json
import urllib2
import urllib
import base64
import gzip
import sys
import requests

__epi__ = "Licence: GPLv3 by Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>"
__version__ = '0.0.1'


API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)


def create_request(request_str, token=API_TOKEN, gzip=False):
    ''' Helper method that formats a URL request for the API by
        wrapping it with your token
    '''
    request = urllib2.Request(request_str)
    base64string = base64.encodestring('%s:%s' % (token, ''))\
        .replace('\n', '')
    request.add_header("Authorization", "Basic %s" % base64string)
    if gzip:
        request.add_header('Accept-encoding', 'gzip')
    return request


def update_BN():

    ''' Main method that downloads EnteroBase data and prepares to load
            into BN.
        '''
    global args
    logging.info('Reading config from %s...' % args.config_file.name)
    CONFIG = get_config(args, args.config_file)
    logging.info('Creating Temporary Config file...')
    format_config(CONFIG)
    for BN_database in CONFIG['databases']:
        if args.fetch_data:
            entero_db = CONFIG['databases'][BN_database]
            output_dir = os.path.join(CONFIG['datadir'], BN_database)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            logging.info('Updating database %s ...' % BN_database)

            logging.info('Fetching %s scheme info from EnteroBase...'
                         % entero_db)
            scheme_dict = get_scheme_info(output_dir, CONFIG['server'],
                                          entero_db, BN_database)

            logging.info('Fetching %s scheme info from EnteroBase...'
                         % entero_db)
            field_list = get_strain_fields(output_dir, CONFIG['server'],
                                           entero_db, BN_database,
                                           scheme_dict.keys())

            logging.info('Fetching %s strain info from EnteroBase...'
                         % entero_db)
            get_strain_data(output_dir, CONFIG['server'], entero_db,
                            BN_database, scheme_dict.keys(), field_list, daily_dump_only = args.daily_dump_only)

#        run_daemon(CONFIG['bnloc'], CONFIG['homedir'], BN_database, CONFIG['script'])

def run_daemon(bnloc, homedir, database, script):

    cmd =  [bnloc, 'homedir=%s' % homedir, 'database=%s' % database, 'script=%s' % script]
    bn_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    bn_proc.wait()
    #logfile = open(os.path.join(os.path.dirname(script), 'error_log.txt'), 'w')
    #logfile.close()

    #logfile = open(os.path.join(os.path.dirname(script), 'error_log.txt'))
    #for log_line in log_printer(logfile):
        #print log_line.strip()
        #if log_line.strip() == 'Done':
            #sys.exit(0)
import StringIO

def get_strain_data(data_dir, server_address , database_name, BN_database, scheme_list, field_list, daily_dump_only=False):
    ''' Fetches strain and ST data.

    '''
    strain_file =  os.path.join(data_dir, '%s-strain_file.gz' %BN_database)
    address = server_address + '/download_straindata?species=%s' % database_name
    response = requests.get(address)
    # Check if any strain has been updated.
    strain_rec = {}
    line_count = 0
    #if os.path.exists(strain_file):
    compressedFile = StringIO.StringIO()
    compressedFile.write(response.content)
    compressedFile.seek(0)

    decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')

    #with  gzip.open(response.raw, 'rb' ) as strain_out:
    for line in decompressedFile.readlines():
        if len(line) > 0 :
            try:
                strain = json.loads(line)
                line_count += 1
                if strain.get('strain_barcode'):
                    if strain_rec.has_key(strain['strain_barcode']):
                        logging.error('Strain barcode already  exists for %s' %strain)
                    strain_rec[strain['strain_barcode']] = strain
                else:
                    logging.error('Strain barcode does not exists for %s' %strain)
            except:
                logging.exception('Could not load line %s' %line)
    logging.info('Read  %d lines' %line_count )
    update_list = []
    logging.info('Loaded %d strain from file' %len(strain_rec.keys()))
    with  gzip.open(os.path.join(data_dir, '%s-strain_file.gz' %BN_database), 'w' ) as strain_out:
        if not daily_dump_only:
            address = server_address + '/api/v2.0/%s/straindata?limit=1000&assembly_status=Assembled,legacy&reldate=1&substrain=true&only_fields=strain_barcode' %database_name
            while True:
                start_time = time.time()
                logging.info('Connecting to %s' %address)
                response = urllib2.urlopen(create_request(address))
                data = json.load(response)
                for strain in data['straindata'].values():
                    strain_barcode =  strain.get('strain_barcode')
                    if strain_barcode is not None:
                        update_list.append(strain_barcode)
                    else:
                        logging.error('No strain barcode for %s' %strain)
                address = data['links']['paging'].get('next')
                if not address:
                    break
                logging.info('Time elasped: %f' % float(time.time() - start_time))
            logging.info('Removed %d outdated strains' %len(update_list) )
            logging.info('Unique %d  strains' %len(set(update_list)) )
            update_list = list(set(update_list))
            found_strains = []
            logging.info('Wrote %d  strains to file ' %len(strain_rec.keys()) )
        # URL for EnteroBase strain records. Fetch 1 only.
            chunks = [update_list[x:x+100] for x in xrange(0, len(update_list), 100)]
            offset = 0
            for chunk in chunks:
                start_time = time.time()
                address = server_address + '/api/v2.0/%s/straindata?assembly_status=Assembled,legacy&substrain=true&barcode=%s&limit=100000' % (database_name, ','.join(chunk))
                logging.info('Connecting to %s' %address)
                response = urllib2.urlopen(create_request(address))
                # Strain data is returned as json.
                data = json.load(response)
                logging.info('Reading %d / %d' %(offset, len(update_list)))
                offset +=   len(chunk)
                seen = {}
                for strain in data['straindata'].values():
                    new_strain = {}
                    if seen.get(strain.get('strain_barcode')):
                        logging.error('Duplicate entry for  %s' %strain.get('strain_barcode'))
                    else:
                        seen[strain.get('strain_barcode')] = True
                    valid_fields = set(field_list) & set(strain.keys())
                    for key in valid_fields:
                        if strain.get(key):
                            new_strain[key] = strain.get(key)
                    if not strain.get('sts'):
                        logging.error('No ST information for %s' %strain.get('strain_barcode'))

                    new_strain.update(_clean_up_st(database_name, strain))

                    if new_strain.get('strain_barcode'):
                        new_strain['lastmodified'] = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%dT%H:%M:%S.%f')
                        strain_rec[new_strain['strain_barcode']] = new_strain
                        found_strains.append(new_strain['strain_barcode'])
                    else:
                        logging.error('Strain has no strain barcode %s ' % new_strain)
                logging.info('Time elasped: %f' % float(time.time() - start_time))
            logging.info('Following strains could not be found: %s' %( ','.join(set(update_list) - set(found_strains))))
        for new_strain in strain_rec.values():
            if new_strain.has_key('sts'):
                new_strain.update(_clean_up_st(database_name, new_strain))
                new_strain.pop('sts', None)
            strain_out.write('%s\r\n' %json.dumps(new_strain))
        strain_out.flush()
        strain_out.close()

def _clean_up_st(database_name, strain):
    new_strain = { }
    for sts in strain.get('sts'):
        if sts.get('st_complex'):
            new_strain['%s_eBG' %sts['scheme_name']] =  sts.get('st_complex')
        if sts.get('lineage'):
            new_strain['%s_lineage' %sts['scheme_name']] =  sts.get('lineage')
        if sts.get('st_id'):
            new_strain['%s_ST' %sts['scheme_name']] =  sts.get('st_id')
        if sts['scheme_name'] == 'rMLST' and database_name == 'senterica':
            if isinstance(sts.get('predicted_serotype'), str):
                new_strain['predicted_serotype'] = sts.get('predicted_serotype')
            else:
                logging.info('Invalid predicted_serotype for %s in %s'
                             %(sts.get('scheme_name'), strain['strain_name']))
    return new_strain

def _get_date(string_time):
    return datetime.datetime.strptime(string_time.split('+')[0], '%Y-%m-%dT%H:%M:%S.%f')


def get_config(args, config_file):
    CONFIG = json.loads(config_file.read())
    required_keys = ['script', 'server', 'bnloc', 'homedir', 'datadir']
    for key in required_keys:
        if getattr(args, key):
            CONFIG[key] = getattr(args, key)
    return CONFIG

def get_scheme_info(data_dir, server_address, database_name, BN_database):
    ''' Creates fields and experimental data types,
        using one strain record as a template.
        Returns list of scheme names, required later.
    '''
    # output files
    scheme_out = open(os.path.join(data_dir, '%s-scheme_file.txt' %BN_database), 'w' )

    # URL for EnteroBase strain records. Fetch 1 only.
    address = server_address + '/api/v2.0/%s/schemes' %database_name
    logging.info('Connecting to %s' %address)
    response = urllib2.urlopen(create_request(address))

    # Strain data is returned as json.
    data = json.load(response)
    scheme_list = []
    for scheme in data['Schemes']:
        if 'MLST' in scheme['scheme_name'] and scheme['scheme_name'] != 'rMLST':
            logging.info('Downloading ST profiles for %s...' %scheme['scheme_name'])
            scheme_list.append(scheme['scheme_name'] )

            response = requests.get(scheme['download_sts_link'], stream=True)
            st_file =  os.path.join(data_dir, '%s-%s-ST.gz' %(BN_database, scheme['scheme_name']))
            with open(st_file, 'wb') as st_file_handle:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        st_file_handle.write(chunk)
                        st_file_handle.flush()
    scheme_dict = {}
    for scheme in scheme_list:
        params = dict(limit=50000, only_fields = 'locus')
        address = server_address + '/api/v2.0/%s/%s/loci?%s' % (database_name, scheme, urllib.urlencode(params))
        logging.info('Connecting to %s' %address)
        response = urllib2.urlopen(create_request(address))
        data = json.load(response)
        scheme_dict[scheme] = []
        for loci in data['loci']:
            scheme_dict[scheme].append(loci['locus'])
    scheme_out.write(json.dumps(scheme_dict))
    return scheme_dict


def get_strain_fields(output_dir, server_address, database_name, BN_database, scheme_list):
    ''' Fetches fields and experimental data types,
        Writes to file to Sync with BN.
    '''
    # output files
    field_names = open(os.path.join(output_dir, '%s-fields_file.txt' %BN_database), 'w' )

    # URL for EnteroBase strain records. Fetch 1 only.
    address = server_address + '/api/v2.0/%s/strains?limit=1' %database_name
    logging.info('Connecting to %s' %address)
    response = urllib2.urlopen(create_request(address))

    # Strain data is returned as json.
    data = json.load(response)
    record = data["Strains"][0]
    field_list = record.keys()
    for scheme in scheme_list:
        field_list.append('%s_ST' %scheme.replace('.','_').replace(' ','_'))
        field_list.append('%s_eBG' %scheme.replace('.','_').replace(' ','_'))
    field_list.append('Predicted_serovar')
    field_list.append('strain_barcode')
    field_list.append('assembly_status')
    field_names.write(json.dumps(field_list))
    return field_list


def format_config(CONFIG):
    ''' I do not know how to pass parameter to the script Engine in BN.
        Dirty hack where this is saved to a file next to script and read in the
        BN script itself later.
    '''
    config_file = os.path.join(os.path.dirname(CONFIG['script']), 'config.json')
    with open(config_file, 'w') as out_config:
        out_config.write(json.dumps(CONFIG))
    out_config.close()


def log_printer(file_handle):
    import time
    file_handle.seek(0,2) # Go to the end of the file
    while True:
        line = file_handle.readline()
        if not line:
            time.sleep(0.1)
            continue
        yield line

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

        parser = argparse.ArgumentParser(description=desc,
                                         epilog=__epi__)
        parser.add_argument('-v', '--verbose', action='store_true',
                            default=False, help='verbose output')

        parser.add_argument('--bnloc', action='store',
                            help='Location of Bionumerics executable',
                            default=None)
        parser.add_argument('--homedir', action='store',
                            help='Bionumerics home directory',
                            default=None)
        parser.add_argument('--datadir', action='store',
                            help='Directory to save EnteroBase data',
                            default=None)
        parser.add_argument('config_file', action='store',
                            help='Config_file', type=file,
                            default='config_file.json')
        parser.add_argument('--script', action='store',
                            help='Script to run', default=None)
        parser.add_argument('--server', action='store',
                            help='EnteroBase Server address', default=None)
        parser.add_argument('--fetch_data', action='store_true',
                            help='Download data from EnteroBase',
                            default=False)
        parser.add_argument('--daily_dump_only', action='store_true',
                            help='Fetch only daily_dump data from EnteroBase',
                            default=False)
        parser.add_argument('--version', action='version',
                            version='%(prog)s ' + __version__)
        args = parser.parse_args()
        if args.verbose:
            print "Executing @ " + time.asctime()

        update_BN()

        if args.verbose:
            print "Ended @ " + time.asctime()
        if args.verbose:
            print 'total time in minutes:',
        if args.verbose:
            print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)