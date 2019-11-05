""" Dependencies:

    iterools
    prody
    glob
    multiprocessing.Pool
    os
    compress_pickle

"""

from predictor.core import training
from predictor.general import pickle_saver

radius = 4
number_of_processors = 27
path_training_folders = "/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"

def main(radius, number_of_processors, path_training_folders):
    training_data = training.training_with_multiprocessing(radius, number_of_processors, path_training_folders)
    filename = "trained_model/UEP_trained_model"
    pickle_saver.saving_pickle(training_data, filename)

if __name__ == "__main__":
    main(radius, number_of_processors, path_training_folders)
