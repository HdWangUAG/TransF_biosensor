import json
import re
import pandas as pd
import os


# 1.1 analysis log file, and extract the affinity

def analysis_log_file(log_path):

    column_names = ['Mode', 'Affinity (kcal/mol)', 'Intramol (kcal/mol)', 'CNN Pose Score', 'CNN Affinity']
    error_df_row = {column: None for column in column_names}
    
    try:
        with open(log_path, 'r' ) as f:
            data = f.readlines()
        if not data:
            return pd.DataFrame([error_df_row]),f"Empty file: {log_path}"
    except Exception as e:
        return pd.DataFrame([error_df_row]),f"Error reading file: {log_path} : {e}"


    # select the indices to process
    
    #-----------------------------------------------------------------------------------------#
    '''adjust the indices here, for rigid docking, it is 19, for flexible docking, it is 20 '''
    #-----------------------------------------------------------------------------------------#
    indices_to_process = [19] 

    extracted_data = []
    error_message = None

    # iterate over selected indices
    for idx in indices_to_process:
        if idx < len(data) and data[idx].strip():
            values = re.findall(r"[-+]?\d*\.\d+|\d+", data[idx])
            
            values = [float(num) if '.' in num else int(num) for num in values]

            if len(values) == 5: # if the length is 5, then it is the correct format
                extracted_data.append(values)
            else:
                print(f'Error in the data: {values}, expected 5 values, but got {len(values)} at index {idx}')
        else:
            error_message = f"Index {idx} out of range or line is empty in file {log_path}"

    # using a dataframe to store the data
    if not extracted_data:
        return pd.DataFrame([error_df_row]), error_message or f"No data collected from {log_path}"

    return pd.DataFrame(extracted_data, columns=column_names), None


# 1.2 for each folder, process all the log files
def process_folder(folder_Path):
    log_files = os.listdir(folder_Path)
    
    combined_df = []
    error_records = []    

    for log_file in log_files:
        

        file = log_file

        data_df, error_msg = analysis_log_file(os.path.join(folder_Path, file))
        data_df['file'] = file
        combined_df.append(data_df)
       
        if error_msg:  # This file had an issue
            error_records.append({'file': file, 'error': error_msg})
    
    #combine all the dataframes
    all_df = pd.concat(combined_df, ignore_index=True)
    error_df = pd.DataFrame(error_records)

    error_df.to_csv(f'{folder_Path}_errors.csv', index=False)
    all_df.to_csv(f'{folder_Path}.csv', index=False)


def load_data():
    for i in range(0,9+1):
        log_path = f'/home/hdwang/protein_sensor_hd/1flm/1flm_1st_docking/repack_conformer{i}/log.csv'

        df_log = pd.read_csv(log_path)

        df = pd.DataFrame()
        # process each row
        df['affinity'] = df_log['Affinity (kcal/mol)']
        df['gnina_score'] = df_log['CNN Pose Score'] * df_log['CNN Affinity']
        df['file'] = df_log['file'].str.split('_', expand=True)[[0, 1,2]].agg('_'.join, axis=1)


        df.to_csv(f'affinity_gnina_score{i}.csv', index=False)

        print(f'Finished processing {log_path}')

'''

dealing with missing value, fill  NaN with the median of each column

'''
# Function to load data and handle missing values
def load_and_process_data(file_paths):
    dfs = []
    for path in file_paths:
        df = pd.read_csv(path)
        # handle missing values
        numeric_columns = df.select_dtypes(include=['number']).columns
        medians = df[numeric_columns].median()
        df[numeric_columns] = df[numeric_columns].fillna(medians)
        dfs.append(df)
    
    # Concatenate all DataFrames into one
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

# Function to calculate aggregation
def aggregate_data(df):
    # Group by 'file' column to sum scores and average affinity
    grouped = df.groupby('file', as_index=False).agg({
        'gnina_score': 'mean',
        'affinity': 'mean'
    })
    return grouped

# Main function to execute the workflow
def main():
    # Define paths to the input files
    file_indices = range(10) # Number of files
    file_paths = [f'affinity_gnina_score{i}.csv' for i in file_indices]
    
    # Process data
    combined_df = load_and_process_data(file_paths)
    
    # Aggregate data
    processed_df = aggregate_data(combined_df)
    
    # Export the processed data to CSV
    processed_df.to_csv('combined_affinity_gnina_score.csv', index=False)
    print("Data has been processed and saved successfully.")




# Execute the main function
if __name__ == "__main__":
    # load the data, mainly log files
    for i in range(0,9+1):
        folder_path = f'/home/hdwang/protein_sensor_hd/1flm/1flm_1st_docking/repack_conformer{i}/log' # D:\4.Testosterone\1st_redesign_docking\repack_conformer0
        df = process_folder(folder_path)
        print(f'Finished processing {folder_path}')
    load_data()
    main()