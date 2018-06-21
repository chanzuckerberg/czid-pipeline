import unittest

import os
import subprocess
import time

from .idseq_step_setup import IdseqStepSetup
from idseq_pipeline.engine.pipeline_flow import PipelineFlow
from idseq_pipeline.steps.run_cdhitdup import PipelineStepRunCDHitDup

class RunCDHitDupTest(unittest.TestCase):

    def test_step_paired(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunCDHitDup, "cdhitdup_out", paired=True)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder

    def test_step_single(self):
        runstep = IdseqStepSetup.get_step_object(PipelineStepRunCDHitDup, "cdhitdup_out", paired=False)
        runstep.start()
        runstep.wait_until_finished()
        # Check results
        # Clean up the folder



