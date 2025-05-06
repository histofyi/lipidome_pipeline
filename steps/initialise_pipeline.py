from typing import Dict, List, Union, Tuple
import os

from shared import load_config


def initialise_pipeline():

    config, success, errors = load_config()

    if success:
        print("Configuration loaded successfully.")

        for path in config['paths']:
            this_path = config["paths"][path]   
            if not os.path.exists(this_path):
                os.makedirs(this_path)
                print(f"Created directory: {this_path}")
            else:
                print(f"Directory already exists: {this_path}")


if __name__ == "__main__":
    initialise_pipeline()
    print("Pipeline initialised.")


