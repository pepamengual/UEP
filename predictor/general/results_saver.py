def saving_file(skempi_predictions, filename):
    with open(filename, "w") as f:
        for name, ratio in sorted(skempi_predictions.items()):
            saving_line = "{} {}".format(name, ratio)
            f.write(saving_line + "\n")
