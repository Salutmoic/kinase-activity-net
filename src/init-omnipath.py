#!/bin/env python2

import pypath

if __name__ == "__main__":
    pa = pypath.PyPath()
    pa.init_network() #exclude=["psite", "signor", "dbptm", "depod", "phelm"])
    # remove high-throughput data
    pa.remove_htp()
    pa.get_directed()
    pa.save_network()
