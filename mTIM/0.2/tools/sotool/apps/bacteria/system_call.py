"""
helper module for pythongrid
"""

import os


def call(command):
    """
    perform system call
    """
    os.system(command)

    return "done"


def call2(arg_tuple):
    """
    performs two system calls
    """

    command1, command2 = arg_tuple

    print command1
    print command2

    os.system(command1)
    os.system(command2)

    return "done"
