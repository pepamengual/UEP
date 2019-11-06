import argparse
from predictor.general import check_file_exists
from predictor.core import training
from predictor.general import pickle_saver
from predictor.general import pickle_reader
from predictor.core import scoring
from predictor.general import results_saver
from predictor.core import scan_file
from predictor.general import save

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
    parser.add_argument('--model', type=str, help='Path of a UEP pre-trained model', default='trained_model/UEP_trained_model')
    parser.add_argument('--skempi', help='Re-run skempi results', action='store_true')
    parser.add_argument('--scan', type=str, help='Scan the entire interface of a given PDB file. A .csv file will be generated', default="")
    args = parser.parse_args()
    return args.cpu, args.radius, args.model, args.skempi, args.scan

def main(number_of_processors=27, radius=4, path_trained_model="trained_model/UEP_trained_model", skempi=False, scan=""):
    path_training_folders="/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"
    
    if not skempi and scan == "":
        raise ValueError("Please, re-run UEP.py using --skempi or --scan arguments")
    else:
        path_trained_model = "{}_4".format(path_trained_model)
        path_trained_model_exists = check_file_exists.file_checker(path_trained_model)
        if path_trained_model_exists:
            training_data = pickle_reader.reading_pickle(path_trained_model)
        else:
            training_data = training.training_with_multiprocessing(radius, number_of_processors, path_training_folders)
            pickle_saver.saving_pickle(training_data, path_trained_model)
        
        if skempi and scan == "":
            print("Running skempi benchmark using {} cpus...".format(number_of_processors))
            skempi_predictions = scoring.scoring_with_multiprocessing(radius, number_of_processors, training_data)
            uep_results_file = "skempi/uep_predictions.txt"
            results_saver.saving_file(skempi_predictions, uep_results_file)
        
        if not skempi and scan != "":
            ratio_dict = scan_file.scan_interface(scan, training_data)         
            save.save_file(ratio_dict, scan)

if __name__ == "__main__":
    cpu, radius, model, skempi, scan = parse_args()
    main(number_of_processors=cpu, radius=radius, path_trained_model=model, skempi=skempi, scan=scan)
