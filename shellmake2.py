import yaml
import argparse

parser = argparse.ArgumentParser(description="Shell script generator for moving FAA files")
parser.add_argument('--group_yml', required=True, help="Path to groups-strains information yaml file")
parser.add_argument('--save_path', required=True, help="Path to save FAA files")
parser.add_argument('--dir_name', required=False, default="faas/", help="Name of the directory to save FAA files")
parser.add_argument('--script', required=False, default="move_faas.sh", help="Name of the shell script. Default is move_faa.sh")
args = parser.parse_args()

# Specify the path to your YAML file
groups_original = args.group_yml
save_path = args.save_path
dir_name = args.dir_name
script = args.script

# Read the YAML file and load it into a Python dictionary
with open(groups_original, 'r') as stream:
    try:
        groups_info = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Some logs for your help
print()
print("This is the group-strain info you provided:")
print(groups_info)
print()
print("Now making shell scripts for moving FAA files into group folders...")

f = open(script, "w")
new_command = ""
main_path = "/jupyterdem/annotation/"

for group in groups_info.keys():
    group_folder = save_path + group + '/' + dir_name
    f.write("mkdir -p" + " " + group_folder + '\n')
    for strain in groups_info[group]:
        faa_path = main_path + group + '_' + strain + '/' + group + '_' + strain + '.faa'
        new_command = "cp" + " " + faa_path + " " + group_folder
        f.write(new_command + '\n')
f.close()

print()
print("Done")