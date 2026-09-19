"""Reject vacuous comparisons and retain all three result comparisons."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from netCDF4 import Dataset

SCRIPT = Path(__file__).resolve().parents[1] / 'verification/internal/compare_outputs.py'


class CompareOutputsTests(unittest.TestCase):
    def test_comparisons(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old, new = root / 'old', root / 'new'
            def run(expected):
                result = subprocess.run([sys.executable, str(SCRIPT), str(old), str(new), str(root / 'report.json')],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
            run(1)  # Neither directory exists.
            old.mkdir()
            new.mkdir()
            run(1)  # Empty directories.
            for path in (old, new):
                (path / 'state.checkpoint').write_bytes(b'committed history')
                (path / 'history.csv').write_text('time,value\n1,2\n')
            run(1)  # Files exist but no numeric arrays are compared.
            for path in (old, new):
                with Dataset(path / 'results.e', 'w') as data:
                    data.createDimension('nodes', 2)
                    data.createVariable('temperature', 'f8', ('nodes',))[:] = [300, 301]
            run(0)
            (new / 'state.checkpoint').unlink()
            run(1)
            (new / 'state.checkpoint').write_bytes(b'changed history')
            run(1)
            (new / 'state.checkpoint').write_bytes(b'committed history')
            (new / 'history.csv').write_text('time,value\n1,3\n')
            run(1)
            (new / 'history.csv').write_text('time,value\n1,2\n')
            with Dataset(new / 'results.e', 'a') as data:
                data['temperature'][1] = 302
            run(1)
            (new / 'results.e').unlink()
            new.rename(root / 'removed')
            run(1)  # Missing current directory with a populated baseline.


if __name__ == '__main__':
    unittest.main()
