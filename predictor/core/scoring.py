def scoring_with_multiprocessing(radius, number_of_processors, training_data):
    from multiprocessing import Pool
    from predictor.core import environment_scoring
    import glob

    skempi_predictions = {}
    
    path_skempi_models = "skempi/foldx_wildtype_models/WT_*"
    pool = Pool(processes=number_of_processors)
    multiple_results = []

    for skempi_pdb_file in glob.glob(path_skempi_models):
        multiple_results.append(pool.apply_async(environment_scoring.score_skempi, (training_data, skempi_pdb_file, radius)))
    for result in multiple_results:
        name, ratio = result.get()
        skempi_predictions.setdefault(name, ratio)
    pool.terminate()

    return skempi_predictions
