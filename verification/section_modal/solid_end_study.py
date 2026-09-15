"""Four literal production cards, full-domain comparison, and reproducible artifact hashes."""
import argparse
import csv
import hashlib
from pathlib import Path

from compare import compare, summary, failed_metrics
from run import run


def study(executable, directory, reference_root):
    source = Path(__file__).resolve().parent
    records = []
    for case in ('axial', 'bending', 'transverse', 'nonuniform'):
        target = directory / case
        name = case + '_hybrid_12'
        mesh = 'plate_solid_ends_128_32.e'
        reference = reference_root / case / (case + '_support_solid.e')
        run(executable, source, target, name, [], mesh)
        measured = compare(target, case + '_hybrid', [12], mesh,
                           reference_path=reference, diagnostics=True, point_support=True)[0]
        failures = failed_metrics(measured)
        report = summary(target / (name + '_summary.csv'))
        records.append(dict(**measured, accuracy_status='FAIL' if failures else 'PASS',
                            failed_metrics=';'.join(failures),
                            global_system_size=int(report['global_system_size']),
                            local_system_size=int(report['local_system_size']),
                            equilibrium_relative_residual=report['equilibrium_relative_residual'],
                            executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                            input_sha256=hashlib.sha256((source / (name + '.fsi')).read_bytes()).hexdigest(),
                            mesh_sha256=hashlib.sha256((source / mesh).read_bytes()).hexdigest(),
                            reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest()))
        with (directory / 'study.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
    if any(row['accuracy_status'] != 'PASS' for row in records):
        raise SystemExit('One or more complete-domain 1% comparisons failed; see study.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('reference_root', type=Path)
    args = parser.parse_args()
    study(args.executable.resolve(), args.directory.resolve(), args.reference_root.resolve())
