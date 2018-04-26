#!/usr/local/env python3

import collections
import os.path
import glob
import json
import time

__author__ = 'duceppemo'


class CheckPhasterServer(object):

    def __init__(self, args):
        # Define variables based on supplied arguments
        self.args = args
        self.folder = args.folder

        # Shared data structure(s)
        self.jobs_dict = collections.defaultdict(dict)

        # run the script
        self.run()

    def run(self):
        self.parse_json()
        self.update_json()
        self.get_ranks()
        self.choose_next()

    def parse_json(self):
        """
        Parse json file into a dictionary
        :return: An updated dictionary
        """
        # Look for all json files in folder
        for json_file in glob.iglob(self.folder + '/*.json'):
            name = os.path.basename(json_file).split('.')[0].split('_')[0]  # Everything before the 1st underscore

            # Parse json file
            data = json.load(open(json_file))
            job_id = data['job_id']
            status = data['status']

            self.jobs_dict[name]['job_id'] = job_id
            self.jobs_dict[name]['status'] = status

            if status in 'Complete':
                zip_url = data['zip']
                self.jobs_dict[name]['zip_url'] = zip_url

    def update_json(self):
        """
        Update all the json files in a folder if results not already downloaded
        :return:
        """

        # Check if some result files are already downloaded
        zip_list = [os.path.basename(z).split('.')[0].split('_')[0] for z in glob.iglob(self.folder + '/*.zip')]

        # Only update the json files with no zip file (Phaster analysis completed and results downloaded)
        for sample, item in self.jobs_dict.items():
            if sample not in zip_list:
                url = 'http://phaster.ca/phaster_api?acc=' + item['job_id']
                print("Updating %s.json" % sample)
                print(url)
                self.download_file(url, self.folder, sample + '.json')

        # Parse the newly downloaded json files
        self.parse_json()

    def get_ranks(self):
        """
        Check if samples analysis on Phaster server is completed or not.
        If not, check how many submissions are ahead of yours

        {"job_id":"ZZ_b7c387e704","status":"5 submissions ahead of yours..."} -> rank = 5
        {"job_id":"ZZ_f86d1fb67c","status":"You're next!..."} -> rank = 2
        {"job_id":"ZZ_f86d1fb67c","status":"Running..."} -> rank = 1
        {"job_id":"ZZ_097f311c05","status":"Complete","zip":"phaster.ca/submissions/ZZ_097f311c05.zip" -> rank = 0
        """

        for sample_name in self.jobs_dict:
            rank = self.jobs_dict[sample_name]['status'].split(" ")[0]
            if rank.isdigit():
                self.jobs_dict[sample_name]['rank'] = int(rank)
            elif rank in 'next':
                self.jobs_dict[sample_name]['rank'] = 2
            elif rank in 'Complete':
                self.jobs_dict[sample_name]['rank'] = 0
            else:  # running
                self.jobs_dict[sample_name]['rank'] = 1

    def choose_next(self):
        """
        Choose if you should download and/or wait
        It takes about 3 min per genome to finish on the Phaster server

        Loops back to itself
        """

        rank_list = [item['rank'] for k, item in self.jobs_dict.items()]
        min_rank = min(rank_list)
        max_rank = max(rank_list)

        completed_samples = [k for k, item in self.jobs_dict.items() if item['rank'] == 0]
        # No need to retrieve that ones we already have
        completed_samples_to_get = [s for s in completed_samples
                                    if not os.path.exists(self.folder + '/' + s + '.zip')]

        # Go ahead and download whatever is completed if don't have it already
        if completed_samples_to_get:
            for sample in completed_samples_to_get:
                url = 'http://' + self.jobs_dict[sample]['zip_url']
                print("Downloading Phaster results for %s" % sample)
                print(url)
                self.download_file(url, self.folder, sample + '_phaster.zip')  # TODO -> Do in a separate thread
        else:  # Wait
            wait_time = max_rank * 3
            print("There are %d submission(s) ahead your fist one." + "\n" +
                  "There are %d submission(s) ahead your last one." + "\n" +
                  "Waiting %d minutes (based on your last submission) before trying to download again..."
                  % (min_rank, max_rank, wait_time))
            time.sleep(wait_time * 60)
            self.update_json()
            self.get_ranks()
            self.choose_next()

    def download_file(self, url, path, name):
        """
        Download a file from an url.
        :param url: url for a specific job_id on Phaster server
        :param path: Where to save the file on disk
        :param name: Name of the file on disk
        :return:
        """

        import urllib.request
        import shutil

        try:
            opener = urllib.request.build_opener()
            opener.addheaders = [('User-agent', 'Mozilla/5.0')]
            opener.open(url)
            # req = urllib.request.Request(url)
            # req.add_header('User-Agent', 'Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11')
            # with urllib.request.urlopen(url) as response, open(path + '/' + name, 'wb') as out_file:
            with opener.open(url) as response, open(path + '/' + name, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        except urllib.error.HTTPError as e:
            print("Error %d. Could not download %s from \"%s\"" % (e.code, name, url))

        time.sleep(1)  # Phaster server seems to be happier by waiting a bit between requests

if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Check status of submissions on Phaster\'s server'
                                        'and download results when available')
    parser.add_argument('-f', '--folder', metavar='/phaster/',
                        required=True,
                        help='Folder holding the json files')

    # Get the arguments into an object
    arguments = parser.parse_args()

    CheckPhasterServer(arguments)


