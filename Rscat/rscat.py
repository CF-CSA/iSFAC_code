"""
 read X-ray scattering factors from a text file, 
 Cromer-Mann parameters and information about the element type
 (Z0, delta Z) and compute Rscat
"""
import argparse
import numpy as np

"""
read scattering factors from a file. Should only contain numbers 
first column is s, other columns for elements
"""
def readfs(filename):
    with open(filename, "r") as f:
        s_fs_data = np.loadtxt(filename)
    return s_fs_data

"""
compute f(s) by the Mott-Bethe formula
Z0, deltaZ are scalars. s, sf are columns
"""
def MB(Z0, deltaZ, s, sf):
    mb = 0.023934 * (Z0 + deltaZ - sf) / np.square(s)
    return mb

def Rscat(fmb, fcm):
    rscat = np.abs(fmb-fcm)/np.abs(fmb)
    return rscat

"""
return Fcalc from CM parametrisation as
f(s) = sum 0..3 a_i*exp(-b_i * s**2) + c
"""


def Fcalc(cm, s):
    fs = (
        cm[8]
        + cm[0] * np.exp(-cm[1] * s**2)
        + cm[2] * np.exp(-cm[3] * s**2)
        + cm[4] * np.exp(-cm[5] * s**2)
        + cm[6] * np.exp(-cm[7] * s**2)
    )
    return fs


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    parser = parser = argparse.ArgumentParser(
        description="compute rscat from list of scattering factors and Cromaer Mann parameters"
    )
    parser.add_argument(
        "-s", "--CM_params", help="Cromer Mann parameters from SFAC card", type=str
    )
    parser.add_argument(
        "-f", "--filename", help="two column file with values for s and for f_X(s)"
    )
    parser.add_argument("-Z", "--Z0", help="value of Z0 for element", type=int)
    parser.add_argument(
        "-z",
        "--deltaZ",
        help="value of Delta Z for element( e.g. C-: -1, Si(IV): +4)",
        type=int,
    )
    parser.add_argument("-S", "--ranges", help="string with lower and upper limit of s", type=str)
    args = parser.parse_args()

    cm = [float(f) for f in (args.CM_params).split()]
    data = readfs(args.filename)
    s_f= data[:, [0, args.Z0]]

    fmb = MB(args.Z0, args.deltaZ, s_f[:, 0], s_f[:,1])
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
