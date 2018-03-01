#!/usr/bin/env python
'''Test of NI-DM SPM export tool

To run the test, copy the NIDM turtle file 'spm_nidm.ttl' obtained by exporting
the results of example001 (as specified in 'examples/spm/example001') in a
directory named spmexport' under 'test' and from the command line call:

python test/TestSPMResultDataModel.py

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>, Satrajit Ghosh
@copyright: University of Warwick 2014
'''
import unittest
import os
from rdflib.graph import Graph
import logging
import glob

# Save debug info in a log file (debug.log)
logging.basicConfig(filename='debug.log', level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)
logger.info(' ---------- Debug log ----------')

from nidmresults.test.test_results_doc import TestResultDataModel
from nidmresults.test.test_commons import *
from nidmresults.test.check_consistency import *
from ddt import ddt, data

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# Find all test examples to be compared with ground truth
test_files = glob.glob(os.path.join(TEST_DIR, 'spmexport', 'ex*', '*.ttl'))
# For test name readability remove path to test file
# test_files = [x.replace(TEST_DIR, "") for x in test_files]
logging.info("Test files:\n\t" + "\n\t".join(test_files))


@ddt
class TestSPMResultsDataModel(unittest.TestCase, TestResultDataModel):

    def setUp(self):
        self.my_execption = ""

        # Display log messages in console
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s')

        gt_dir = os.path.join(TEST_DIR, 'spmexport', 'ground_truth')

        TestResultDataModel.setUp(self, gt_dir)

    # @data(*test_files)
    # def test_class_consistency_with_owl(self, ttl):
    #     """
    #     Test: Check that the classes used in the ttl file are defined in the
    #     owl file.
    #     """
    #     ex = self.load_graph(ttl)
    #     ex.owl.check_class_names(ex.graph, ex.name, True)

    # @data(*test_files)
    # def test_attributes_consistency_with_owl(self, ttl):
    #     """
    #     Test: Check that the attributes used in the ttl file comply with
    # their
    #     definition (range, domain) specified in the owl file.
    #     """
    #     ex = self.load_graph(ttl)
    #     ex.owl.check_attributes(ex.graph, ex.name, True)

    @data(*test_files)
    def test_examples_match_ground_truth(self, ttl):
        """
        Test03: Comparing that the ttl file generated by SPM and the expected
        ttl file (generated manually) are identical
        """
        ex = self.load_graph(ttl)
        ex_gt = ""

        # Creating a single ground truth graph (by merging all included ground
        # truths)
        gt = Graph()
        run_with_octave = False
        first = True

        if 'full_example' in ttl:
            gt_oct = Graph()
        else:
            gt_oct = None

        for gt_file in ex.gt_ttl_files:
            logging.info("Ground truth ttl: " + gt_file)

            # RDF obtained by the ground truth export
            gt.parse(gt_file, format='turtle')

        # TODO: Workaround to deal with differences in shasum for gzipped files
        # between octave and Matlab and across hosts
        # -- Ignore shasums --
        CRYPT_RDFLIB = rdflib.Namespace(
            "http://id.loc.gov/vocabulary/preservation/cryptographicHashFunctions#")
        gt.remove((None, CRYPT_RDFLIB.sha512, None))
        ex.exact_comparison = True
        # End of TODO

        self.compare_full_graphs(gt, ex.graph, ex.owl,
                                 ex.exact_comparison, False)

        ex_gt += self.my_execption

        if ex_gt:
            raise Exception(ex_gt)


if __name__ == '__main__':
    unittest.main()
