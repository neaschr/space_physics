import pandas as pd

def data_reader(datafile_name: str, column_names: list) -> pd.DataFrame:
    '''Function that takes in the name a data file and the names of the columns and returns
    the data as a pandas dataframe'''

    #Reading the data from the file given 
    data = pd.read_csv(datafile_name, comment = '%', header = None, delim_whitespace = True, names = column_names)

    return data

