""" Dependencies:

    iterools
    prody
    glob
    multiprocessing.Pool
    os
    compress_pickle

"""
from predictor.general import check_file_exists
from predictor.core import training
from predictor.general import pickle_saver
from predictor.general import pickle_reader
from predictor.core import scoring
from predictor.general import results_saver

radius = 4
number_of_processors = 27
path_training_folders = "/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"

def main(radius, number_of_processors, path_training_folders):
    path_trained_model = "trained_model/UEP_trained_model"
    path_trained_model_exists = check_file_exists.file_checker(path_trained_model)
    if path_trained_model_exists:
        training_data = pickle_reader.reading_pickle(path_trained_model)
    else:
        training_data = training.training_with_multiprocessing(radius, number_of_processors, path_training_folders)
        pickle_saver.saving_pickle(training_data, path_trained_model)
    
    path_skempi_models = "/home/pepamengual/UEP/skempi/foldx_wildtype_models/WT_*"
    skempi_predictions = scoring.scoring_with_multiprocessing(radius, number_of_processors, training_data, path_skempi_models)
    
    uep_results_file = "skempi/uep_predictions.txt"
    results_saver.saving_file(skempi_predictions, uep_results_file)

if __name__ == "__main__":
    main(radius, number_of_processors, path_training_folders)
