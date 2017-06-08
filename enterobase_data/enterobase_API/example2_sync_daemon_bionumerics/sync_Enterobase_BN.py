#! /usr/bin/env python
"""
Syncronises BN database with EnteroBase cache in local directory

### CHANGE LOG ###
2016-08-30 Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>
    * Initial build

"""

import os
import subprocess
import time
import traceback
import logging
import bns
import urllib2, json, time, base64
import gzip
import md5

__epi__ = "Licence: GPLv3 by Nabil-Fareed Alikhan <n-f.alikhan@warwick.ac.uk>"
__version__ = '0.0.1'


def load_params(config_file):

	with open(config_file) as handle:
		CONFIG = json.loads(handle.read())
		return CONFIG

def update_field_names(field_file):
	if 'is_uberstrain' not in bns.Database.Db.Fields:
		bns.Database.Db.Fields.Add('is_uberstrain')
	with open(field_file) as field_info:
		for line in field_info.readlines():
			field_names = json.loads(line)
			field_names.remove('created')
			field_names.remove('uberstrain')
			for field in field_names:
				if field not in bns.Database.Db	.Fields:
					logging.info('Added %s ' %field)
					bns.Database.Db.Fields.Add(field)
		field_names.append('rMLST_eBG')
		field_names.append('rMLST_ST')
		bns.Database.Db.Fields.Save()
		return field_names

def update_schemes(scheme_file):
	scheme_list  = []
	with open(scheme_file) as scheme_info:
		for line in scheme_info.readlines():
			schemes = json.loads(line)
			for scheme in schemes:
				scheme_list.append(scheme)
				try :
					bns.Database.ExperimentType(scheme)
				except :
					bns.Database.Db.ExperimentTypes.Create(scheme,'CHR')
					logging.info( 'Adding new fields for %s' %scheme)
					curr_scheme = bns.Characters.CharSetType(scheme)
					if scheme != 'wgMLST':
						for loci in schemes[scheme]:
							if curr_scheme.FindChar(loci) == -1:
								curr_scheme.AddChar(loci, 100000)
	return scheme_list

def update_exp_data(scheme, scheme_file):
	st_update = {}
	checksum_field  = '%s_checksum' %scheme
	loci = []
	exp = bns.Characters.CharSet()
	experiment = bns.Characters.CharSetType(scheme)
	for id in xrange(experiment.GetCount()):
		loci.append(experiment.GetChar(id))
	# Create a checksum field if it does not exist
	if checksum_field not in bns.Database.Db	.Fields:
		bns.Database.Db.Fields.Add(checksum_field)
	# Build a lookup of all strains of the same STs (for bulk insertion)
	for strain_barcode in get_selected_keys():
		strain_barcode = str(strain_barcode)
		if not strain_barcode or strain_barcode == 'None':
			continue
		ST_id = bns.Database.EntryField(strain_barcode, '%s_ST' %scheme ).Content
		if st_update.get(ST_id):
			st_update[ST_id].append(strain_barcode)
		else:
			st_update[ST_id] = [strain_barcode]
	logging.info('Loaded %d profiles' %len(st_update))
	# Read gz file for each scheme
	if scheme_file.endswith('.gz'):
		headers = gzip.open(scheme_file).readlines()[0].strip().split('\t')
		file_handle = gzip.open(scheme_file).readlines()[1:]
	else:
		headers = open(scheme_file).readlines()[0].strip().split('\t')
		file_handle = open(scheme_file).readlines()[1:]
	count = 0
	for line in file_handle:
		count += 1
		if count % 1000 == 0 :
			logging.info('Updated %d STs' %count)
		allele_profile = dict(zip(headers, line.strip().split('\t')))
		ST = allele_profile['ST']
		new_checksum = _create_checksum(allele_profile)
		# Check if allele profile has been modified via checksums
		# If so update
		if  not st_update.get(ST, None):
			logging.info('No record with ST %s in scheme : %s ' %(ST, scheme))
		first_update = True
		for strain_barcode in st_update.get(ST, [] ):
			if bns.Database.EntryField(strain_barcode, '%s_eBG' %scheme ).Content == 'None':
				bns.Database.EntryField(strain_barcode, '%s_eBG' %scheme ).Content = ''
			if int(ST) < 1:
				bns.Database.EntryField(strain_barcode, '%s_checksum' %scheme ).Content  = ''
				bns.Database.EntryField(strain_barcode, '%s_ST' %scheme ).Content = 'Invalid ST'
				bns.Database.EntryField(strain_barcode, '%s_eBG' %scheme ).Content = ''
				bns.Database.Experiment(strain_barcode, scheme).Delete()
			else:
				old_checksum = bns.Database.EntryField(strain_barcode, '%s_checksum' %scheme ).Content
				if (new_checksum !=  old_checksum) or old_checksum == '' or old_checksum is None:
					if first_update:
						data = [float(allele_profile.get(fld, 0)) for fld in loci]
						presence = [ int(d>0) for d in data ]
						first_update = False
					if not exp.Load(strain_barcode, scheme) :
						exp.Create(scheme, '', strain_barcode)
					else:
						d2, p2 = [], []
						exp.GetData(d2, p2)
						if not cmp(d2, data) :
							continue
					bns.Database.EntryField(strain_barcode, '%s_checksum' %scheme ).Content  = _create_checksum(allele_profile)
					exp.SetData(data, presence)
					exp.Save()
		bns.Database.Db.Fields.Save()

def _create_checksum(allele_profile):
	allele_profile.pop('ST', None)
	new_checksum = [value for (key, value) in  sorted(allele_profile.items())]
	checksum_md5 = md5.new(','.join(new_checksum))
	return checksum_md5.hexdigest()

def update_strain_meta(strain_file, field_names):
	keys = []
	for ent in bns.Database.Db.Entries:
		keys.append(str(ent))
	logging.info('added %d keys' %len(keys))
	with gzip.open(strain_file) as strain_info:
		found = []
		good_read = False
		for line in strain_info.readlines():
			try:
				if len(line) > 0 :
					strains_metadata = json.loads(line)
					found.append(str(strains_metadata['strain_barcode']))
			except Exception:
				good_read = False
				logging.exception(line)
				logging.exception('Error reading line, not valid json: %s' %line)
		if good_read:
			no_result = __diff(keys, found)
			if len(no_result) > 0:
				bns.Database.Db.Entries.Delete(no_result)
		found = list(set(found))
		logging.info('added %d strains' %len(found))
		add_strains = __diff(found, keys)
		logging.info('Number of new records %d strains' %len(add_strains))
		if len(add_strains) > 0:
			bns.Database.Db.Entries.Add(add_strains)
		logging.info('Added Strains')
	with gzip.open(strain_file) as strain_info:
		found = []
		count = 0
		for line in strain_info.readlines():
			try:
				strain = json.loads(line)
				count += 1
				strain_barcode = str(strain.get('strain_barcode'))
				if strain_barcode == strain['uberstrain']:
					bns.Database.EntryField(strain_barcode, 'is_uberstrain').Content = 'True'
				else:
					bns.Database.EntryField(strain_barcode, 'is_uberstrain').Content = 'False'
					logging.info('Setting %s to substrain' %strain_barcode)
				for field in field_names:
					try:
						if strain.get(field) and strain.get('strain_barcode'):
			#				val = str(strain.get(field).replace(u"\uFFFD","").replace(u"\ufffd","")).encode('ascii', 'ignore')
							val = str(strain.get(field))
							if val != 'None' and val is not  None:
								
								bns.Database.EntryField(str(strain.get('strain_barcode')), field).Content = str(val)
					except:
						logging.error('Error loading field %s in %s' %(field, strain_barcode))
				if count % 1000 == 0 or count == len(found):
					logging.info('Updated %d records' %count)
					bns.Database.Db.Fields.Save()							
			except Exception:
				logging.error('Could not encode %s ' %strain_barcode)			
		bns.Database.Db.Fields.Save()

def get_selected_keys():
	keys = []
	all_keys = [] 
	for ent in bns.Database.Db.Entries:
		if ent.Selection == 1:
			keys.append(str(ent))		
		all_keys.append(str(ent))
	if len(keys) > 0 : 
		return keys
	else: 
		return all_keys

def __diff(a, b):
	b = set(b)
	return [aa for aa in a if aa not in b and aa is not None]


DEBUG = True
EXP_ONLY  = False
if __name__ == '__main__' :
	if DEBUG:
		script_dir = os.path.dirname(os.path.realpath(__file__))
	else:
		script_dir = os.path.dirname(os.path.realpath("__file__"))
	logging.basicConfig(filename=os.path.join(script_dir, 'error_log.txt'),level=logging.INFO)
	CONFIG = load_params(os.path.join(script_dir, 'config.json'))
	data_file_prefix = os.path.join(CONFIG['datadir'], bns.Database.Db.Info.Name)
	field_names = update_field_names(os.path.join(data_file_prefix, '%s-fields_file.txt' %bns.Database.Db.Info.Name))
	scheme_list = update_schemes(os.path.join(data_file_prefix, '%s-scheme_file.txt' %bns.Database.Db.Info.Name))
	new_field = [] 
	for field in field_names : 
		if field.endswith('_ST') or field.endswith('eBG'):
			new_field.append(field)
	if EXP_ONLY: 
		field_names = new_field
	update_strain_meta(os.path.join(data_file_prefix, '%s-strain_file.gz' %bns.Database.Db.Info.Name), field_names)
	scheme_list = ['MLST_Achtman', 'rMLST', 'cgMLST', 'wgMLST']
	for scheme in scheme_list:
		update_exp_data(scheme, os.path.join(data_file_prefix, '%s-%s-ST.gz' %(bns.Database.Db.Info.Name, scheme)) )
	logging.info('Done')
