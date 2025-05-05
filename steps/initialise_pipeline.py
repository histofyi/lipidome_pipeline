from typing import Dict, List, Union, Tuple
import os

from shared import load_json, load_config


def initialise_pipeline():

    config, success, errors = load_config()

    if success:
        print("Configuration loaded successfully.")

        for path in config["paths"]:
            if not os.path.exists(path):
                os.makedirs(path)
                print(f"Created directory: {path}")
            else:
                print(f"Directory already exists: {path}")


if __name__ == "__main__":
    initialise_pipeline()
    print("Pipeline initialised.")


