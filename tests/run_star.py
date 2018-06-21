import unittest

from .idseq_step_setup import IdseqStepSetup
from idseq_dag.steps.run_star import PipelineStepRunStar

import pdb


class RunStarTest(unittest.TestCase):
    def test_step_paired(self):
        run_step = IdseqStepSetup.get_step_object(PipelineStepRunStar, "star_out", paired=True)
        pdb.set_trace()
        run_step.start()
        run_step.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder

    def test_step_single(self):
        run_step = IdseqStepSetup.get_step_object(PipelineStepRunStar, "star_out", paired=False)
        pdb.set_trace()
        run_step.start()
        run_step.wait_until_finished()
        # TODO: Check results
        # TODO: Clean up the folder
