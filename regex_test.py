#!/usr/bin/env python
# coding=utf-8

import re
import argparse


def main():
    parser = argparse.ArgumentParser(description="for help with regex, visit: https://regex101.com/r/F561kR/3")

    parser.add_argument('-t', '--text', type=str,
                        default="TR34869|C1_G2_I2|M.23583_R0",
                        help='test text to search')

    parser.add_argument('-r', '--regex', type=str, default=r'^(.+?)\|(.+?)_(.+?)_.+$',
                        help='test pattern to search with')

    args = parser.parse_args()

    print("\ncreating key from:\n{}\n\nwith:\n{}\n".format(args.text, args.regex))

    search_obj = re.search(args.regex, args.text)
    key = search_obj.groups()

    print("key:\n{}\n\n".format(key))


if __name__ == "__main__":
    main()
