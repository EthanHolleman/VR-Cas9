import json

# Read json file path
def read_config_json(config_path):
    with open(str(config_path)) as handle:
        config  = json.load(handle)
    return config


def adjust_paths(config_dict, path_dict):
    pass


def rewrite_config():
    pass


def main():
    
    pass

if __name__ == '__main__':
    main()