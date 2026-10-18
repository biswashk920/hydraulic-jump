"""Runs the whole validation suite (gvf/validate.py) and requires every check to pass."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gvf import validate


class TestValidationSuite(unittest.TestCase):
    def test_every_check_passes(self):
        items = validate.run_all()
        self.assertGreater(len(items), 30)
        failed = [f"{i['group']}: {i['name']} (measured {i['measured']:.2e}, limit {i['limit']:.0e})" for i in items if not i["ok"]]
        self.assertEqual(failed, [])


if __name__ == "__main__":
    unittest.main()
