import os
import sys
from argparse import ArgumentParser

from bancotetrado.structure import Structure
from bancotetrado import banco_painter



def DrawFromFile(filename_json, output_file):
    structure = Structure().fromFile(filename_json)
    return banco_painter.Draw(structure, output_file)




def main():
    
    parser = ArgumentParser("bancotetrado")
    parser.add_argument("-i", "--input", required = True, help = "path to JSON file")
    parser.add_argument("-o", "--output", help = "output path")

    args = parser.parse_args()

    root, _ = os.path.splitext(args.input)
    output_template = args.output
    if output_template == None:
        output_template = root

    DrawFromFile(args.input, output_template)


if __name__ == "__main__":
    sys.exit(main())
