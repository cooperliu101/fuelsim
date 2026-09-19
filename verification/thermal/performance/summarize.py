"""Summarize measured samples without rerunning either solver."""
import argparse
import json
from pathlib import Path
import statistics


def summarize(work):
    measurements = json.loads((work/'timing_summary.json').read_text())
    result = {}
    for case,data in measurements.items():
        accuracy = json.loads((work/'validation'/case/'comparison.json').read_text())
        assert accuracy['passed'] and len(data['samples']) == 3
        entry = {'accuracy':accuracy, 'fuelsim_external_seconds':data['fuelsim_external'],
                 'abaqus_external_seconds':data['abaqus_external'],
                 'abaqus_over_fuelsim_external_ratio':data['abaqus_over_fuelsim_external_ratio']}
        for name,solver,key in [('fuelsim_workflow_seconds','fuelsim','total_seconds'),
                                ('abaqus_analysis_wall_seconds','abaqus','analysis_wall_seconds'),
                                ('abaqus_analysis_cpu_seconds','abaqus','analysis_cpu_seconds')]:
            values = [sample[solver][key] for sample in data['samples']]
            entry[name] = {'median':statistics.median(values),'min':min(values),'max':max(values)}
        entry['fuelsim_jacobian_evaluations'] = data['samples'][0]['fuelsim']['jacobian_evaluations_total']
        entry['abaqus_factorizations'] = data['samples'][0]['abaqus']['factorizations']
        entry['maximum_pointwise_relative_error_percent'] = 100*max(
            field.get('maximum_pointwise_relative',0) for field in accuracy['fields'].values())
        result[case] = entry
    (work/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('work',type=Path)
    summarize(parser.parse_args().work)
