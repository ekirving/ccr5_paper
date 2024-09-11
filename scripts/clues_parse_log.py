#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json

import click
from scipy.stats.distributions import chi2


@click.command()
@click.option("--rsid", metavar="<string>", help="RefSeq ID", required=True)
@click.option("--ancestry", metavar="<string>", help="Ancestral path", required=True)
@click.option("--use-freq", metavar="<string>", help="Was the modern frequency used", required=True)
@click.option("--mod-freq", metavar="<string>", help="The modern population frequency", required=True)
@click.option("--log", "log_file", metavar="<file>", type=click.Path(writable=True), help="Log file", required=True)
@click.option("--out", "output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def clues_parse_log(rsid, ancestry, use_freq, mod_freq, log_file, output):
    """
    Parse the Clues log file to extract the information we want.
    """
    with open(log_file) as fin:
        epochs = dict()

        for line in fin:
            if "logLR" in line:
                lnl_ratio = line.split().pop()
            elif "epoch" in line and "selection" in line:
                # handle multiple epochs
                while True:
                    try:
                        epoch, s = next(fin).split()
                        epochs[epoch] = float(s)
                    except ValueError:
                        break

        # convert the log-likelihood ratio into a p-value
        # https://en.wikipedia.org/wiki/Wilks%27_theorem
        pval = chi2.sf(2 * float(lnl_ratio), 1)

        data = {
            "rsid": rsid,
            "mode": "ancient",
            "ancestry": ancestry,
            "use_freq": use_freq,
            "mod_freq": mod_freq,
            "logLR": float(lnl_ratio),
            "pval": pval,
            "epochs": epochs,
        }

    json.dump(data, output, indent=2)


if __name__ == "__main__":
    clues_parse_log()
