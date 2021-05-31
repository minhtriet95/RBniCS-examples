# A class to automatically clean old offline training data
#
# Created: 2020
# Last modify: May 31, 2021

import os
import shutil


class Cleaner():
    def __init__(self, foldername=str):
        self.foldername = foldername
        self.path = os.getcwd()

    # Get all the folder names in current working directory
    def get_names(self, path):
        return os.listdir(self.path)

    def delete_offline_folder(self):
        items = self.get_names(self.path)
        # Check if the offline folder is exist
        if self.foldername in items:
            shutil.rmtree(self.path + "/" + self.foldername)
        else:
            raise NameError(f"Cannot find folder name: {self.foldername}")
        