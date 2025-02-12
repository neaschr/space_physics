import pandas as pd

def data_reader(datafile_name: str) -> pd.DataFrame:

    #Defining the column names for the data
    column_names = [
        'Height_km',
        'Oksygen-atoms_cm-3',
        'Nitrogen-molecules_cm-3',
        'Oksygen-molecules_cm-3',
        'Mass_density_g_per_cm3',
        'Temperature_neutral_K'
    ]

    #Reading the data from the file given 
    data = pd.read_csv(datafile_name, comment = '%', header = None, delim_whitespace = True, names = column_names)

    return data