import unittest

from idseq_pipeline.engine.pipeline_flow import PipelineFlow
import pdb


class CustomTest1(unittest.TestCase):
    def test_custom1(self):
        # TODO: Check results
        # TODO: Clean up the folder

        flow = PipelineFlow(lazy_run=False,
                            dag_json="examples/custom_test_1.json")
        flow.start()
