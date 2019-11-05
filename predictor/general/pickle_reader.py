def reading_pickle(path_to_file):
    from compress_pickle import load
    training_data = load(path_to_file, compression="lzma", set_default_extension=False)
    return training_data
