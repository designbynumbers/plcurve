import unittest
import os

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

suite = unittest.TestLoader().discover(os.path.dirname(__file__))
if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite)
