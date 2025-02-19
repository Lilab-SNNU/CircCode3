# -*- coding = utf-8 -*-
from .utils.argument_parser import init_argparse


def main():
    parser = init_argparse()
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()