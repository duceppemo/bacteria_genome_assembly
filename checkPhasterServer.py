#!/usr/local/env python3

import collections
import os.path
import glob
import json
import time
import os
import sys


__author__ = 'duceppemo'
__version__ = '0.2.3'


class CheckPhasterServer(object):

    def __init__(self, args):
        # Define variables based on supplied arguments
        self.args = args
        self.input_folder = args.input
        self.output_folder = args.output
        self.check = args.check
        self.submit = args.submit

        # Shared data structure(s)
        self.fasta_list = list()
        self.fasta_dict = collections.defaultdict(dict)
        self.jobs_dict = collections.defaultdict(dict)

        # run the script
        self.run()

    def run(self):
        # check which flag has been used
        if self.submit:
            if not self.input_folder:
                self.error('Please provide an input folder with fasta files to submit using the "-i" option')
            # Work on fasta files
            self.list_assemblies(self.input_folder)

            # Setup progress bar
            counter = 0
            total_samples = len(self.fasta_list)
            for fasta in self.fasta_list:
                if self.check_fasta(fasta) is True:
                    counter += 1
                    sample_name = os.path.basename(fasta).split('.')[0].split('_')[0]
                    sys.stdout.write('\rSubmitting sample "' + sample_name
                                     + '" (' + str(counter) + '/' + str(total_samples) + ')')
                    self.submit_assembly(fasta)
                else:
                    print("The following file does now meet PHASTER requirements for minimum contig length" +
                          " and wont be submitted for phage analysis:")
                    print(fasta + "\n")
                # Clear dictionary
                self.fasta_dict.clear()
                sys.stdout.flush()
            print('\n')

        if self.check:
            if not self.output_folder:
                self.error('Please provide an output folder with json files to check status using the "-o" option')
            # Work on json files retrieved from the phaster server
            self.parse_json()
            self.update_json()
            self.get_ranks()
            self.choose_next()

        if not self.check and not self.submit:
            self.error("Please use the %s and/or the %s flag(s)" %('"--check"', '"--submit"'))

    def parse_fasta(self, f, d):
        """
        Parse fasta file into dictionary
        :param f: Input fasta file
        :param d: empty dictionary
        :return: populated dictionary
        """

        seqid, seq = None, []
        with open(f, 'r') as file:
            for line in file:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seqid and seq:
                        d[seqid] = ''.join(seq)
                        seqid = line[1:].split(" ")[0]
                        seq = []
                    else:  # for the first entry
                        seqid = line[1:].split(" ")[0]
                else:
                    seq.append(line)
            if seqid and seq:  # for the last entry
                d[seqid] = ''.join(seq)

    def list_assemblies(self, folder):
        """
        List the fasta files to upload to phaster server.
        :param folder: The folder containing the assemblies in fasta format.
        :return: A populated list of fasta file
        """

        for assembly in glob.iglob(folder + '/*.fasta'):
            self.fasta_list.append(assembly)

        if not self.fasta_list:
            self.error('The provided input folder contains no ".fasta" files')

    def check_fasta(self, f):
        """
        Check if fasta files meets phaster server requirements
        if one entry:
            min_contig_size = 1500 bp
        elif 2 entries or more:
            min_contig_size = 2000 bp
        else:
            error(no entry in fasta file
        :param f: input fasta file
        :return: True or False
        """

        self.parse_fasta(f, self.fasta_dict)
        flag = None

        if len(self.fasta_dict.keys()) == 1:
            if len(self.fasta_dict.values()) >= 1500:
                flag = True
            else:
                flag = False
        elif len(self.fasta_dict.keys()) == 0:
            flag = False
        else:  # > 1
            for id, seq in self.fasta_dict.items():
                if len(seq) >= 2000:
                    flag = True
                else:
                    flag = False
                    break
        return flag

    def submit_assembly(self, f):
        """
        Submit the fasta file to the server via their API
        :return:
        """

        import requests

        if len(self.fasta_dict.keys()) == 1:  # if one contig
            api_url = 'http://phaster.ca/phaster_api'
        else:  # if more than 1 contig
            api_url = 'http://phaster.ca/phaster_api?contigs=1'

        sample = os. path.basename(f).split('_')[0]
        with open(f, 'rb') as my_file:
            data = my_file.read()

        # Header must have 'Content-Type': 'application/x-www-form-urlencoded', thus need to use "data:
        # "files" wont work, "multipart/form-data" is not supported
        r = requests.post(api_url, data=data)  # response is a json file
        # prepared = r.request  # Debug

        # if status not OK (a 4XX client error or 5XX server error response)
        # if r.status_code != 200:  # 200 = OK
        if r.status_code != requests.codes.ok:
            # Write error to file
            with open(self.output_folder + '/' + sample + '.error', 'w') as json_out_handle:
                json_out_handle.write(r.text)
            # r.raise_for_status()  # Not sure I Want to raise and Exception here...
        else:
            # save json to file
            with open(self.output_folder + '/' + sample + '.json', 'w') as json_out_handle:
                json_out_handle.write(r.text)

    def parse_json(self):
        """
        Parse json file into a dictionary
        :return: An updated dictionary
        """

        # Look for all json files in folder
        for json_file in glob.iglob(self.output_folder + '/*.json'):
            name = os.path.basename(json_file).split('.')[0].split('_')[0]  # Everything before the 1st underscore

            # Parse json file
            try:
                data = json.load(open(json_file))
                job_id = data['job_id']
                self.jobs_dict[name]['job_id'] = job_id
                try:
                    status = data['status']
                    self.jobs_dict[name]['status'] = status
                    if status in 'Complete':
                        zip_url = data['zip']
                        self.jobs_dict[name]['zip_url'] = zip_url
                except KeyError:
                    self.jobs_dict[name]['job_id'] = job_id
                    self.jobs_dict[name]['status'] = 'error'
                    print("There was a problem running PHASTER for %s (%s)" % (name, job_id))
            except json.decoder.JSONDecodeError:
                print("JSON file for sample %s is empty. Upload process to PHASTER server failed." % name)

    def update_json(self):
        """
        Update all the json files in a folder if results not already downloaded
        :return:
        """

        from time import sleep

        # Check if some result files are already downloaded
        zip_list = [os.path.basename(z).split('.')[0].split('_')[0] for z in glob.iglob(self.output_folder + '/*.zip')]

        counter = 0
        to_update = len(self.jobs_dict.keys()) - len(zip_list)
        # Only update the json files with no zip file (Phaster analysis completed and results downloaded)
        if to_update > 0:
            print("\n%d json files to update:" % to_update)
            for sample, item in self.jobs_dict.items():
                if sample not in zip_list:
                    counter += 1
                    url = 'http://phaster.ca/phaster_api?acc=' + item['job_id']
                    print("Updating %s.json (%d/%d)" % (sample, counter, to_update))
                    print(url)
                    self.download_file(url, self.output_folder, sample + '.json')
                    sleep(10)
        else:
            print("All result (zip) files already downloaded, exiting")
            sys.exit()  # exit script, nothing else to do

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
        {"job_id":"ZZ_ccc1cf07f5","error":"There was a problem running PHASTER..."} ->
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
                                    if not os.path.exists(self.output_folder + '/' + s + '_phaster.zip')]

        # Go ahead and download whatever is completed if don't have it already
        if completed_samples_to_get:
            print("\n%d new samples ready to download:" % len(completed_samples_to_get))
            for sample in completed_samples_to_get:
                url = 'http://' + self.jobs_dict[sample]['zip_url']
                print("Downloading Phaster results for %s" % sample)
                print("\t%s" % url)
                self.download_file(url, self.output_folder, sample + '_phaster.zip')  # TODO -> Do in a separate thread
        if max_rank > 0:  # Wait
            wait_time = max_rank * 3
            remaining = len(list(filter(lambda x: x > 0, rank_list)))  # how many ranks not equals to zero
            print("\nYou still have %d samples waiting for processing:"
                  % remaining)
            print("\tThere are %d submission(s) ahead of you before your first submitted file is processed."
                  % min_rank)
            print("\tThere are %d submission(s) ahead of you before your last submitted file is processed"
                  % max_rank)
            print("Waiting %d minutes (based on longer estimated wait time) before trying to download again..."
                  % wait_time)
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
            with opener.open(url) as response, open(path + '/' + name, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        except urllib.error.HTTPError as e:
            print("Error %d. Could not download %s from \"%s\"" % (e.code, name, url))

        time.sleep(10)  # Phaster server seems to be happier by waiting a bit between requests (Authors recommend 60s)

    def error(self, message):
        """
        Display error type and help
        :param message: Error string
        :return: Error and help message
        """
        sys.stderr.write('error: %s\n' % message)
        parser.print_help()
        sys.exit(2)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Submit fasta and/or check status of submissions on Phaster\'s server'
                                        'and download results when available. Will run util all the samples'
                                        'have been retrieved. Status is updated every "n" minutes,'
                                        'where n = number of samples in queue * 3 minutes')
    parser.add_argument('--submit', action='store_true',
                        help='Use this flag to submit your fasta file(s) to phaster'
                             'Requires "-i"'
                             'A json status file will be created after submission')
    parser.add_argument('--check', action='store_true',
                        help='Check status of submission on phaster server.'
                             'Requires "-o".'
                             'Looking for json files in the output folder ("-o")'
                             'Will update json files and download result file if analysis completed')
    parser.add_argument('-i', '--input', metavar='/assembly_folder/',
                        required=False,
                        help='Folder holding the input fasta files')
    parser.add_argument('-o', '--output', metavar='/result_folder/',
                        required=True,
                        help='Folder to save the phaster result files'
                             'Required'
                             'Also used to save the json status files from the phaster API')

    # Get the arguments into an object
    arguments = parser.parse_args()

    CheckPhasterServer(arguments)
