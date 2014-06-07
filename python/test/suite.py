import unittest
import os

suite = unittest.TestLoader().discover(os.path.dirname(__file__))
if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite)
