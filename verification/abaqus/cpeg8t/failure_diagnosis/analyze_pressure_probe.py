"""Report native activation separately from the existing affine conductance law."""
import json
from pathlib import Path
from create_pressure_probe import GAPS

ROOT = Path(__file__).resolve().parent


def analyze():
    steps = json.loads((ROOT / 'pressure_opening_fields.json').read_text())['steps']
    rows = []
    for index, gap in enumerate(GAPS):
        fields = steps[f'GAP_{index}'][-1]['fields']
        reactions = fields['RFL11']['values']
        upper = sum(v['data'] for v in reactions if 9 <= v['nodeLabel'] <= 16)
        lower = sum(v['data'] for v in reactions if 1 <= v['nodeLabel'] <= 8)
        pressure = next(v['values'] for k, v in fields.items() if k.startswith('CPRESS'))
        rows.append(dict(gap_m=gap, native_upper_heat_W=upper, native_lower_heat_W=lower,
                         balance_W=upper+lower,
                         native_maximum_pressure_Pa=max(v['data'] for v in pressure),
                         affine_law_heat_W=.002*100*(1100+.002*max(-1e9*gap, 0))))
    (ROOT / 'pressure_activation.json').write_text(json.dumps(rows, indent=2)+'\n')
    return rows


if __name__ == '__main__':
    print(json.dumps(analyze(), indent=2))
