import runpy
import sys
import os
import filecmp
import shutil

root_disk = os.getcwd()

# Clean before test
for item in os.listdir(root_disk):
    if item not in ['expected_output', 'integrity_test.py', 'PDBs_test']:
        if(os.path.isdir(item)):
            shutil.rmtree(item)
        else:
            os.remove(item)

# Set Machaon's console arguments
os.makedirs('output', exist_ok=True)
sys.argv = ['', '1', root_disk, 'PDBs_test', os.path.join(root_disk, 'output'), '6VXX.A', '1273', '1', '--noThirdPartyData']

# Call Machaon
working_directory = root_disk.replace('/test', '/src')
os.chdir(working_directory)
runpy.run_path(path_name=os.path.join(working_directory, 'run.py'), run_name="__main__")

# Check output
filecmp.dircmp(os.path.join(root_disk, 'output'), os.path.join(root_disk, 'expected_output')).report_full_closure()