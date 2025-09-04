"""
This Python program is an wrapper for executing a command-line tool called "haf-pipe".

@author: Luis Javier Madrigal-Roca

@date: 2025-09-04

"""

import subprocess
import os
import sys
import argparse

def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run haf-pipe with specified parameters.")

    parser.add_argument('--haf_wrapper', type=str, required=True, help='Path to the bash haf-pipe wrapper.')
    parser.add_argument('--harp_exe', type=str, required=True, help='Path to the harp executable.')

    

if __name__ == "__main__":
    main()