"""The  accessionid2seq command."""
from .base import Base
import sys
import time
import subprocess
import random
import pickle
import shelve
import argparse
import gzip
import re
import threading
import traceback

class Accessionid2seq(Base):
    def run(self):
        print("Placeholder")
