#!/usr/local/env python3

__author__ = 'duceppemo'
__version__ = '0.1'


import os
from nested_dict import nested_dict


class BDA(object):

    def __init__(self, args):
        import multiprocessing

        self.args = args
        self.input = args.input
        self.out = args.output

        # Number of cpu
        self.cpus = int(multiprocessing.cpu_count())

        self.blast_dict = nested_dict()

        # Default output path is input assembly path
        if not self.out:
            self.out = os.path.dirname(self.input) + "/"

        self.position_list = list()  # contain tuples (contig, position, sense

        # run the script
        self.run()

    def run(self):
        self.check_dependencies()
        self.parse_fasta(self.input)
        self.run_blastn('nt', self.input)

    def parse_fasta(self, f):
        """
        Parse input fasta file using SeqIO from Biopython
        :param f: Assembly file in fasta format
        :return: Populated dictionary
        """
        from Bio import SeqIO

        self.contigs_dict = SeqIO.to_dict(SeqIO.parse(f, format='fasta'))

    def run_blastn(self, ref_db, query):
        """
        Perform blastn using biopython
        :param ref: A fasta file for which "makeblastdb' has already been run
        :param query: The assembly, in fasta format
        :return: XML file
        """
        from Bio.Blast.Applications import NcbiblastnCommandline
        from io import StringIO

        blastn = NcbiblastnCommandline(db=ref_db, query=query, evalue='1e-10',
                                             outfmt=5, max_target_seqs=5,
                                             num_threads=self.cpus)
        (stdout, stderr) = blastn()

        if stderr:
            raise Exception('There was a problem with the blast')

        # Search stdout for matches - if the term Hsp appears (the .find function will NOT
        # return -1), a match has been found, and stdout is written to file
        if stdout.find('Hsp') != -1:
            blast_handle = StringIO(stdout)  # Convert string to IO object
            self.parse_blast_output(blast_handle)
        else:
            print("No reference sequences found.")

        # TODO -> Check that no contigs have both the dnaA and the repA genes
        # Raise a warning to the user if is the case
        
    def identify_MIH(hit_list, bit_scores):
        """
        :param hit_list: a list containing the descriptors of hits for the current qresult
        :param bit_scores: a list of bit scores associated with the hits from the qresult
        :return: the most informative hit (MIH) descriptor
        """
        from nltk.corpus import stopwords
        from sklearn.feature_extraction.text import CountVectorizer
        from sklearn.cluster import KMeans
        import numpy as np

        # Get the informative part of the hit description first by removing stopwords (common words)
        filtered_list = []
        for hit in hits_list:
            words_list = hit.lower().split()
            filtered_words = [word for word in words_list if word not in stopwords.words('english')]
            descriptor = " ".join(filtered_words).split("os=", 1)[0]
            filtered_list.append(descriptor)
        # Convert the description into a vector which can be used for clustering
        vectorizer = CountVectorizer()
        vectorized_hits = vectorizer.fit_transform(filtered_list)
        vector_array = np.array(vectorized_hits.toarray())
        # Create cluster centers, then test descriptions for which cluster they fall into
        if len(hits_list) < 8:
            kmeans = KMeans(n_clusters=len(hits_list), n_init=100).fit(vector_array)
        else:
            kmeans = KMeans(n_init=100).fit(vector_array)
        cluster_associations = kmeans.predict(vector_array)
        # Associate the bit scores with the cluster predictions
        cluster_aggregates = [(cluster, score) for cluster, score in zip(cluster_associations, bit_scores)]
        cluster_scores = [0] * 8
        # Determine which clustering contains the highest score as a combination of number of bit scores and their values
        for cluster, score in cluster_aggregates:
            cluster_scores[cluster] += score
        # Find out what descriptions the highest scored cluster is associated with, and return the description 
        MIH_index = list(cluster_associations).index(cluster_scores.index(max(cluster_scores)))
        return hits_list[MIH_index]

    def parse_blast_output(self, blast_handle):
        """
        test
        :param blast_handle: An xml Blast output file handle from io.StringIO
        :return:
        """
        from Bio import SearchIO
        from Bio.Seq import Seq
        from os import remove

        records_dict = SearchIO.to_dict(SearchIO.parse(blast_handle, 'blast-xml'))
        fixed_seq = None
        
        file = open(self.input, "r+")
        file_text = file.read()
        file.close()
            
        
        # Create the list of words to filter out uninformative hits
        filter_strings = ['putative', 'protein like', 'protein related', 'contains similarity to',
                          'predicted', 'hypothetical protein', 'unnamed protein product',
                          'unknown', 'expressed protein', 'uncharacterized', 'probable',
                          'possible', 'potential']

        for id, qresult in records_dict.items():  
            hit_descriptions = []
            e_values = []
            bit_scores = []
            for hit in qresult.hits:
                hit_desc = hit._description
                # Filter out uninformative hits here, create a list contianing the remaining hits
                if not any([filter in str(hit_desc).lower() for filter in filter_strings]):
                    hit_descriptions.append(hit_desc)
                    e_values.append(hit.hsps[0].evalue)
                    bit_scores.append(hit.hsps[0].bitscore)
            # Assuming there are hits left after the filtering, find the MIH
            # If no hits are left after filtering, just return "hypothetical protein"
            if len(hit_descriptions) > 0:
                proper_desc = identify_MIH(hit_descriptions, bit_scores)
                file_text_array = file_text.split("/n")
                for line in file_text_array:
                    if id in line:
                        file_text_array[file_text_array.index(line)] = ">{1} {2}".format(id, proper_desc)
        file_path = os.path.abspath(self.input)
        os.remove(file_path)
        new_file = open(file_path, "w"):
            new_file.write("\n".join(file_text_array))

        # Close file StringIO handle
        blast_handle.close()

        return fixed_seq

    def check_dependencies(self):
        pass


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Reorder assembly to have the dnaA gene first')
    parser.add_argument('-i', '--input', metavar='input.fasta',
                        required=True,
                        help='Input fasta file to blast')
    parser.add_argument('-o', '--output', metavar='output_bda.fasta',
                        required=True,
                        help='Annotated fasta file with useful descriptions')

    # Get the arguments into an object
    arguments = parser.parse_args()

    BDA(arguments)
