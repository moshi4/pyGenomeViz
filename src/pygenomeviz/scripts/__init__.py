import argparse
from typing import List

from pygenomeviz.genbank import Genbank


def get_argparser(prog_name: str) -> argparse.ArgumentParser:
    """Get argument parser

    Parameters
    ----------
    prog_name : str
        Program name to be used

    Returns
    -------
    parser : argparse.ArgumentParser
        Argument parser
    """

    class CustomHelpFormatter(argparse.RawTextHelpFormatter):
        def __init__(self, prog, indent_increment=2, max_help_position=40, width=None):
            super().__init__(prog, indent_increment, max_help_position, width)

    desc = f"pyGenomeViz CLI utility for visualization workflow using {prog_name}"
    epilog = "[*] marker means the default value."
    parser = argparse.ArgumentParser(
        description=desc,
        epilog=epilog,
        add_help=False,
        formatter_class=CustomHelpFormatter,
    )
    return parser


def print_args(args: argparse.Namespace) -> None:
    """Print arguments

    Parameters
    ----------
    args : argparse.Namespace
        Arguments
    """
    header = "\n" + "*" * 20 + " Run Parameters " + "*" * 20
    print(header)
    for k, v in args.__dict__.items():
        if type(v) == list:
            print(f"{k}:")
            for i, f in enumerate(v, 1):
                print(f"  {i:02d}. {f}")
        else:
            print(f"{k}: {v}")
    print("*" * len(header))


def gbk_files2gbk_objects(gbk_files: List[str]) -> List[Genbank]:
    """Convert genbank files to Genbank objects

    Parameters
    ----------
    gbk_files : List[str]
        Genbank files (Target range can be set as follows 'file:100-1000')

    Returns
    -------
    gbk_list : List[Genbank]
        Genbank objects
    """
    gbk_list: List[Genbank] = []
    for gbk_file in gbk_files:
        contents = str(gbk_file).split(":")
        gbk_file = contents[0]
        if len(contents) == 1:
            gbk_list.append(Genbank(gbk_file))
        else:
            ranges = contents[1].split("-")
            min_range, max_range = int(ranges[0]), int(ranges[1])
            gbk_list.append(Genbank(gbk_file, min_range=min_range, max_range=max_range))
    return gbk_list
