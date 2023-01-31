from __future__ import annotations

import argparse

from pygenomeviz.parser import Genbank


def get_argparser(prog_name: str | None = None) -> argparse.ArgumentParser:
    """Get argument parser

    Parameters
    ----------
    prog_name : str | None
        Program name to be used

    Returns
    -------
    parser : argparse.ArgumentParser
        Argument parser
    """

    class CustomHelpFormatter(argparse.RawTextHelpFormatter):
        def __init__(self, prog, indent_increment=2, max_help_position=40, width=None):
            super().__init__(prog, indent_increment, max_help_position, width)

    desc = "pyGenomeViz CLI utility for visualization workflow"
    if prog_name is not None:
        desc += f" using {prog_name}"
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
        if k != "format" and type(v) == list:
            print(f"{k}:")
            for i, f in enumerate(v, 1):
                print(f"  {i:02d}. {f}")
        else:
            print(f"{k}: {v}")
    print("*" * len(header))


def gbk_files2gbk_objects(gbk_files: list[str]) -> list[Genbank]:
    """Convert genbank files to Genbank objects

    Parameters
    ----------
    gbk_files : list[str]
        Genbank files.
        User can optionally specify genome range and reverse complement.
        - Example1. Set 100 - 1000 range `file:100-1000`
        - Example2. Set reverse complement `file::-1`
        - Example3. Set 100 - 1000 range of reverse complement `file:100-1000:-1`

    Returns
    -------
    gbk_list : list[Genbank]
        Genbank objects
    """
    gbk_list: list[Genbank] = []
    for gbk_file in gbk_files:
        contents = str(gbk_file).split(":")
        gbk_file = contents[0]
        if len(contents) == 1:
            gbk_list.append(Genbank(gbk_file))
        elif len(contents) == 2:
            ranges = contents[1].split("-")
            min_range, max_range = int(ranges[0]), int(ranges[1])
            gbk_list.append(Genbank(gbk_file, min_range=min_range, max_range=max_range))
        else:
            reverse = True if int(contents[2]) == -1 else False
            if contents[1] == "":
                gbk_list.append(Genbank(gbk_file, reverse=reverse))
            else:
                ranges = contents[1].split("-")
                min_range, max_range = int(ranges[0]), int(ranges[1])
                gbk = Genbank(
                    gbk_file, min_range=min_range, max_range=max_range, reverse=reverse
                )
                gbk_list.append(gbk)
    return gbk_list
