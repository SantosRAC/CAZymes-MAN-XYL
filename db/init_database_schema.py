from functions import *
import argparse

parser= argparse.ArgumentParser(description='Download genomes from NCBI')
parser.add_argument('--password', metavar='password', type=str,
                    help='password for the database',
                    required=True)

args= parser.parse_args()

createDB(password=args.password)