import pandas as pd

def readhealthkitcsv(csvpath):

    # Read the CSV file using Pandas, skipping the first 2 lines (MATLAB index starts at 1, hence lines 1-2 are skipped)
    # and parse 'startDate' and 'endDate' columns as datetime, considering timezone
    df = pd.read_csv(
        csvpath,
        skiprows=1,
        sep=',',
    )

    # Parse 'startDate' and 'endDate' using pd.to_datetime() with the format
    df['startDate'] = pd.to_datetime(df['startDate'], format='%Y-%m-%d %H:%M:%S %z', errors='coerce')
    df['endDate'] = pd.to_datetime(df['endDate'], format='%Y-%m-%d %H:%M:%S %z', errors='coerce')

    return df

