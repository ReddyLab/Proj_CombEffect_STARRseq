### get arguments
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--arg1', help='first argument')
parser.add_argument('--arg2', help='second argument')
args=parser.parse_args()
print("arguements")
for key, val in vars(args).items():
    print(f"{key}: {val}")
