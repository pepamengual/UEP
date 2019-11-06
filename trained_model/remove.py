

def remove_training_data():
    from compress_pickle import load

    training_data_path = "/home/pepamengual/UEP/trained_model/UEP_trained_model_4"
    skempi_data_path = "/home/pepamengual/UEP/trained_model/substracted_4"
    substracted_model = {}
    training_data = load(training_data_path, compression="lzma", set_default_extension=False)
    skempi_data = load(skempi_data_path, compression="lzma", set_default_extension=False)

    for environment, amino_acid_dict in training_data.items():
        for amino_acid, counts in amino_acid_dict.items():
            if environment in skempi_data and amino_acid in skempi_data[environment]:
                substract = counts - skempi_data[environment][amino_acid]
                substracted_model.setdefault(environment, {}).setdefault(amino_acid, substract)
            else:
                substracted_model.setdefault(environment, {}).setdefault(amino_acid, counts)
    from compress_pickle import dump

    dump(substracted_model, "substracted_def_4", compression="lzma")

def main():
    remove_training_data()
main()
