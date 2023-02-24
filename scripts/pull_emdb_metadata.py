import pickle
import numpy as np
from ftplib import FTP
from operator import itemgetter

import os
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import coloredlogs, logging
from time import process_time


def get_emdb_entries():
    servername = "files.rcsb.org"
    ftp = FTP(servername)
    ftp.login(user='random_user',passwd='random_password')

    subdir = "pub/emdb/structures"

    try:
        entries_raw = ftp.nlst(subdir)

    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            print("No files in this directory")
        else:
            raise

    entries_cleaned = [x.split('/')[-1] for x in entries_raw if '.' not in x.split('/')[-1]]
    emdb_entries = entries_cleaned
    
    return emdb_entries


def extract_metadata(x):
    data = x['crossreferences']['citation_list']['primary_citation']['citation_type']
    return data

def retrieveInfoEMDB(emdb_id, url_prefix='https://www.ebi.ac.uk/emdb/api/entry/publications/'):    
    emdb_url = url_prefix + emdb_id
    
    # set up session to allow retry 3 times
    # backoff_factor will help to apply delays between attempts
    # avoiding ConnectionError due to request quota
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    
    respond = session.get(url=emdb_url)

    if respond.status_code != 200:
        logging.error("request status %d" % (respond.status_code))
        raise DownloadException("Error, retrieveInfoFromEMDB request status %d" % (respond.status_code))
    else:
        data = respond.json()
    
    try:
        metadata = extract_metadata(data)
        return metadata
    except Exception as e:
        logging.error(e)
        return None

def partition_range(N, num_partitions):
    partition_size = N // num_partitions
    partitions = [range(i * partition_size, (i + 1) * partition_size) for i in range(num_partitions - 1)]
    partitions.append(range((num_partitions - 1) * partition_size, N))
    
    return partitions

def get_metadata(emdb_entries, interval_length=100):
    N_entries = len(emdb_entries)
    num_partitions = int(N_entries/interval_length)
    index_intervals = partition_range(N_entries, num_partitions)

    metadata = {}
    for i in range(len(index_intervals)):
        entries_interval = itemgetter(*index_intervals[i])(emdb_entries)
        try:
            annotations = {emdb_id:retrieveInfoEMDB(emdb_id) for emdb_id in entries_interval}
            metadata.update(annotations)

        except Exception as e:
            print(e)
    
    return metadata

if __name__ == "__main__":
    emdb_entries = get_emdb_entries()
    metadata = get_metadata(emdb_entries[:len(emdb_entries)])

    with open('data/metadata_all_23-Feb-2023.pickle', 'wb') as fp:
        pickle.dump(metadata, fp)
