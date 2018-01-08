"""Tests for our main idseq_pipeline CLI module."""


from subprocess import PIPE, Popen as popen
from unittest import TestCase

from idseq_pipeline import __version__ as VERSION


class TestHelp(TestCase):
    def test_returns_usage_information(self):
        output = popen(['idseq_pipeline', '-h'], stdout=PIPE).communicate()[0]
        self.assertTrue('Usage:' in output)

        output = popen(['idseq_pipeline', '--help'], stdout=PIPE).communicate()[0]
        self.assertTrue('Usage:' in output)


class TestVersion(TestCase):
    def test_returns_version_information(self):
        output = popen(['idseq_pipeline', '--version'], stdout=PIPE).communicate()[0]
        self.assertEqual(output.strip(), VERSION)
