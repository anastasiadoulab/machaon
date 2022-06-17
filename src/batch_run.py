import argparse
import sys
sys.path.append('..')
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Machaon v1.0 : Find common folds between proteins.')
    parser.add_argument('cores', type=int, help='The number of cores available for Machaon.')
    
    args = vars(parser.parse_args())

    config_manager = ConfigurationManager()
    configurations = config_manager.parse_config()

    comparison = Machaon()
    comparison.max_cores = args['cores']
    comparison.perform_comparisons(configurations)
