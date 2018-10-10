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

    def parse_blast_output(self, blast_handle):
        """
        test
        :param blast_handle: An xml Blast output file handle from io.StringIO
        :return:
        """
        from Bio import SearchIO
        from Bio.Seq import Seq

        records_dict = SearchIO.to_dict(SearchIO.parse(blast_handle, 'blast-xml'))
        fixed_seq = None

        for id, qresult in records_dict.items():  # this is a sorted dictionare, with best hit first?
            for hit in qresult.hits:
                for hsps in hit.hsps:
                    query_seq = self.contigs_dict[id].seq

                    """
                    https://groups.google.com/forum/#!topic/blast2go/B6g6PfX6YW8
                    
                    Removed words are: putative, protein like,
                    protein related, contains similarity to, similar to predicted
                    
                    Removed descriptions are:  hypothetical protein, unnamed protein
                    product, unknown, expressed protein
                    """

                    filter_words = ['putative', 'hypothetical']

                    # TODO -> put good info in object to be filtered

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
