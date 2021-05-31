import os
import shutil


class Cleaner():
    def __init__(self, foldername=str):
        self.foldername = foldername
        # Get current folder path
        self.path = os.getcwd()

    # List all objects in current directory
    def list_names(self, path):
        return os.listdir(self.path)

    def delete_offline_folder(self):
        items = self.list_names(self.path)
        # Check if the offline folder is exist
        if self.foldername in items:
            shutil.rmtree(self.path + "/" + self.foldername)
        else:
            raise NameError(f"Cannot find folder name: {self.foldername}")


if __name__ == "__main__":
    clean = Cleaner("GaussianEIM")
    # Return to the upper folder
    clean.path = clean.get_folder_path() + "/.."
    clean.delete_offline_folder()
