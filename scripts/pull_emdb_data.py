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


def extract_data(x):
    data = {}
    for i, y in enumerate(x['structure_determination_list']['structure_determination']):
        data[i] = {}
        for j, z in enumerate(y['specimen_preparation_list']['specimen_preparation']):
            data[i][j] = z
    return data

def retrieveInfoEMDB(emdb_id, url_prefix='https://www.ebi.ac.uk/emdb/api/entry/experiment/'):
    #extract_data = lambda x:x['structure_determination_list']['structure_determination'][0]['specimen_preparation_list']['specimen_preparation'][0]
    
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
        exptl_annotations = extract_data(data)
        return exptl_annotations
    except Exception as e:
        logging.error(e)
        return None

def partition_range(N, num_partitions):
    partition_size = N // num_partitions
    partitions = [range(i * partition_size, (i + 1) * partition_size) for i in range(num_partitions - 1)]
    partitions.append(range((num_partitions - 1) * partition_size, N))
    
    return partitions

def get_exptl_annotations(emdb_entries, interval_length=100):
    N_entries = len(emdb_entries)
    num_partitions = int(N_entries/interval_length)
    index_intervals = partition_range(N_entries, num_partitions)

    exptl_annotations = {}
    for i in range(len(index_intervals)):
        entries_interval = itemgetter(*index_intervals[i])(emdb_entries)
        try:
            annotations = {emdb_id:retrieveInfoEMDB(emdb_id) for emdb_id in entries_interval}
            exptl_annotations.update(annotations)

        except Exception as e:
            print(e)
    
    return exptl_annotations

if __name__ == "__main__":
    emdb_entries = get_emdb_entries()
    exptl_annotations = get_exptl_annotations(emdb_entries[:len(emdb_entries)])

    with open('data/exptl_annotations_all_22-Feb-2023.pickle', 'wb') as fp:
        pickle.dump(exptl_annotations, fp)
