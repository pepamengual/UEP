def save_file(ratio_dict, scan):
    import pandas as pd
    df = pd.DataFrame(ratio_dict)
    output_file = scan.split("/")[-1].split(".pdb")[0]
    df.to_csv("{}_UEP.csv".format(output_file))
