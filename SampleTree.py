#!/complgen/bin/anaconda/bin/python2.7
# wdecoster

import fdb
import pygraphviz as pgv
import getpass
import argparse
import sys
import os


def main():
    user, pw, individuals, noshow = initialize()
    for individual in individuals:
        sampledict = getInfo(user, pw, individual)
        drawTree(sampledict, individual, noshow)


def initialize():
    '''
    Parsing of arguments and getting login information from the user
    '''
    parser = argparse.ArgumentParser(
        description="Create a sample tree with all samples and availability for an individual.")
    parser.add_argument(
        "sample", help="Sample(s) for which a sampletree should be made.", nargs='+')
    parser.add_argument("--noshow", help="Don't show trees after they are made.",
                        action="store_true")
    args = parser.parse_args()
    user = getpass.getpass("Please provide NBD_SC gentli user account:")
    pw = getpass.getpass("Please provide password for gentli user {}:".format(user))
    if not len(args.sample) > 1 and os.path.isfile(args.sample[0]):
        indiv = [line.strip() for line in open(args.sample[0]).readlines() if not line == ""]
    else:
        indiv = args.sample
    return(user, pw, indiv, args.noshow)


def getInfo(user, pw, id):
    '''
    Doing two calls to gentli to get sample information from two tables.
    Information gathered is structured in a dict containing
     key: sample numberSample for which a sampletree should be made
     value: tuple of sampletype and source sample (parent)
    Number of locations found in sample_location dictates how many vials of a sample is available (only applicable to biobank samples)
    '''
    try:
        con = fdb.connect(dsn='molgenvz.cde.ua.ac.be:/home/firebird/gentli.fdb',
                          user=user, password=pw, role='NBD_SC')
        cur = con.cursor()
    except fdb.fbcore.DatabaseError as ConnError:
        print("\n\nCONNECTION ERROR: There was a problem connecting to the database. Perhaps your password or user were incorrect?")
        print("Please try again. If the error persists, please loudly yell for help and show error message below.\n\n.")
        sys.exit(ConnError)
    cur.execute(
        'select "sample", "type", "source_sample" from "NBD_SC:sample" where "individual" = ?', (id,))
    sampledict = {}
    for sample in cur.fetchall():
        sampledict[sample[0]] = (sample[1], sample[2])
    if len(sampledict.keys()) == 0:
        sys.exit(
            "\n\nINPUT ERROR: Invalid individual {}\nIndividual {} not found in NBD_SC:sample.".format(id, id))
    cur.execute('select "sample" from "NBD_SC:sample_location" where "individual" = ?', (id,))
    samples = [item[0] for item in cur.fetchall()]
    for s in set(samples):
        sampledict[s] = (str(samples.count(s)) + "x " + sampledict[s][0], sampledict[s][1])
    return sampledict


def drawTree(sampledict, individual, noshow):
    '''
    Generate a tree by iterating over the dictionary and adding edges.
    If a sample doesn't have a parent (root or orphan) a KeyError will be generated, and the sample will be added as node
    '''
    G = pgv.AGraph(directed=True)
    for sample in sorted(sampledict.keys()):
        try:
            G.add_edge(
                str(sampledict[sample][1]) + r'\n' + str(sampledict[sampledict[sample][1]][0]),
                str(sample) + r'\n' + str(sampledict[sample][0])
            )
        except KeyError:  # In the case that this sample doesn't have a parent, is root or orphan
            G.add_node(str(sample) + r'\n' + str(sampledict[sample][0]))
    G.layout(prog='dot')
    output = 'SampleTree_' + individual + '.png'
    try:
        G.draw(output)
        print("Finished, created file {}".format(output))
    except IOError:
        sys.exit(
            "\n\nPERMISSION ERROR: Could not create graph.\nDo you have writing permission in the current directory?")
    if not noshow:
        os.system('eog ' + output)


if __name__ == "__main__":
    main()
