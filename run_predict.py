import os 

def run_boltz_batch_prediction_sequentially(input_yaml_folder: str):
    """
    Runs Boltz prediction sequentially for all YAML files within a specified folder.

    Args:
        input_yaml_folder (str): The path to the folder containing batch YAML files.

    """
    # Get a list of all YAML files in the input folder
    yaml_files = [f for f in os.listdir(input_yaml_folder) if f.endswith('.yaml')]

    for yaml_file in yaml_files:
        yaml_path = os.path.join(input_yaml_folder, yaml_file)
        # boltz command in shell
        boltz_command = f"boltz predict {yaml_path}  --use_msa_server --recycling_steps 20 --output_format pdb --diffusion_samples 5"
        
        os.system(boltz_command)
        
if __name__ == "__main__":

    input_folder = "/home/hdwang/TransF_biosensor/yaml"

    run_boltz_batch_prediction_sequentially(input_folder)