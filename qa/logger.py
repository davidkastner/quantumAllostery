"""Format log commands"""

import os
import time


def clean_exit():
    """
    Clean exit from the program.
    """

    sys.exit("Thanks for using quantumAllostery.\nExited successfully.\n")

def nonoption_exit():
    """
    Exit the program when an invalid option is selected.
    """

    sys.exit("Sorry that option is not available.\nPlease select one of the available options.")

def help():
    """
    Print the help menu.
    """

    print("> GitHub: https://github.com/davidkastner/quantumAllostery")
    print("> Documenation: https://quantumallostery.readthedocs.io")
