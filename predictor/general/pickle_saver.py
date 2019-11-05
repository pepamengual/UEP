def saving_pickle(data, filename):
    from compress_pickle import dump
    dump(data, filename, compression="lzma", set_default_extension=False)
