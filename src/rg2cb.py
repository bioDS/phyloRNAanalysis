#!/usr/bin/env python3
"""Add RG tag as a CB tag."""
import argparse
import os
import pysam

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Add RG tag as a CB tag to a SAM/BAM file."""
        )
    parser.add_argument("input", type=str, help="an input SAM/BAM file")
    parser.add_argument("output", type=str, help="an output SAM/BAM file")
    args = parser.parse_args()
    return args

def main():
    """Add RG tag as a CB tag to a SAM/BAM file."""
    args = parse_args()

    imod = mode(args.input)
    omod = mode(args.output)

    with pysam.AlignmentFile(args.input, "r"+imod) as inbam, \
         pysam.AlignmentFile(args.output, "w"+omod, template=inbam) as outbam:
        for read in inbam.fetch(until_eof=True):
            rg2cb(read)
            outbam.write(read)

def rg2cb(read):
    """Add RG tag as a CB tag to a read"""
    text = read.get_tag("RG")
    read.set_tag("CB", text, value_type="Z")

def mode(file):
    """Detect pysam SAM/BAM file mode"""
    ext = os.path.splitext(file)[1]
    if ext == ".bam":
        return "b"
    if ext == ".sam":
        return ""
    raise ValueError(f"File {os.path.basename(file)} must be either SAM or BAM")

if __name__ == "__main__":
    main()
