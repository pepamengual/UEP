def file_checker(file_path):
    from os import path
    file_exists = path.isfile(file_path)
    if file_exists:
        print("---> File {} has been found...".format(file_path))
    else:
        print("--> File {} has not been found...".format(file_path))
    return file_exists
