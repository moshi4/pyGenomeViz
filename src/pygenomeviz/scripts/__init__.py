import argparse


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
        def __init__(self, prog, indent_increment=2, max_help_position=30, width=None):
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
    header = "*" * 20 + " Run Parameters " + "*" * 20
    print(header)
    for k, v in args.__dict__.items():
        if type(v) == list:
            print(f"{k}:")
            for i, f in enumerate(v, 1):
                print(f"  {i:02d}. {f}")
        else:
            print(f"{k}: {v}")
    print("*" * len(header))
