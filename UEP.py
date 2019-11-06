""" Dependencies:

    argparse
    itertools
    prody
    glob
    multiprocessing.Pool
    os
    compress_pickle

"""
import argparse
from predictor.general import check_file_exists
from predictor.core import training
from predictor.general import pickle_saver
from predictor.general import pickle_reader
from predictor.core import scoring
from predictor.general import results_saver

number_of_processors = 27
path_training_folders = "/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"

HELP = " \
Command:\n \
----------\n \
run skempi benchmark: python3 UEP.py --cpu 27 --skempi True \n \
scan interface of PDB file: python3 UEP.py --cpu 27 --scan PDB.pdb \
"

def parse_args():
    parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cpu', type=int, help='Amount of CPUs to use along the simulation', default=27)
    parser.add_argument('--radius', type=int, help='PPI interface radius', default=4)
    parser.add_argument('--train', type=str, help='Path to interactome3D training data', default='/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*')
    parser.add_argument('--model', type=str, help='Path of a UEP pre-trained model', default='trained_model/UEP_trained_model')
    parser.add_argument('--skempi', help='Re-run skempi results', action='store_true')
    parser.add_argument('--scan', type=str, help='Scan the entire interface of a given PDB file')
    args = parser.parse_args()
    return args.cpu, args.radius, args.train, args.model, args.skempi, args.scan

def main(number_of_processors=27, radius=4, path_training_folders="/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*", path_trained_model="trained_model/UEP_trained_model", skempi=False, scan=""):
    
    path_trained_model_exists = check_file_exists.file_checker(path_trained_model)
    if path_trained_model_exists:
        training_data = pickle_reader.reading_pickle(path_trained_model)
    else:
        training_data = training.training_with_multiprocessing(radius, number_of_processors, path_training_folders)
        pickle_saver.saving_pickle(training_data, path_trained_model)
    
    skempi_predictions = scoring.scoring_with_multiprocessing(radius, number_of_processors, training_data)
    
    uep_results_file = "skempi/uep_predictions.txt"
    results_saver.saving_file(skempi_predictions, uep_results_file)

if __name__ == "__main__":
    cpu, radius, train, model, skempi, scan = parse_args()
    main(number_of_processors=cpu, radius=radius, path_training_folders=train, path_trained_model=model, skempi=skempi, scan=scan)
