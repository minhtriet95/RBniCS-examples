import os
import shutil


class Cleaner():
    def __init__(self, foldername=str):
        self.foldername = foldername
        # Get current working directory
        self.path = os.getcwd()

    # List all objects in current directory
    def list_all_objs_in_current_dir(self):
        return os.listdir(self.path)

    def delete_folder(self):
        items = self.list_all_objs_in_current_dir()
        # Check if the offline folder is exist
        if self.foldername in items:
            shutil.rmtree(self.path + "/" + self.foldername)
        else:
            raise NameError(f"Cannot find folder name: {self.foldername}")


if __name__ == "__main__":
    clean = Cleaner("UnitCell")  # Remember to change "foldername"
    # Return to the upper folder
    clean.path = clean.path + "/.."
    clean.delete_folder()
