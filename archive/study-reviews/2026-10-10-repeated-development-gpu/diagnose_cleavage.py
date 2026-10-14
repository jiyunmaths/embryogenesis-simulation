"""Descriptive matched-lineage clock/axis differences; no alignment or solver."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path('output/studies/repeated-development-gpu-qualification-v2')
HERE=Path(__file__).parent


def main():
    rows=[];pins={}
    for branch in ('off','tension_adhesion','polarity','full'):
        lineages=[];events=[]
        for level in ('coarse','fine'):
            folder=ROOT/(branch+'_'+level)
            for name in ('gpu-lineage.json','gpu-events.json'):
                path=folder/name;pins[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
            lineages.append(json.loads((folder/'gpu-lineage.json').read_text()))
            events.append({e['parent']:e for e in json.loads((folder/'gpu-events.json').read_text())})
        assert [(r['id'],r['parent']) for r in lineages[0]]==[(r['id'],r['parent']) for r in lineages[1]]
        assert events[0].keys()==events[1].keys()
        records=[]
        for parent,a in events[0].items():
            b=events[1][parent]
            assert a['children']==b['children'] and a['draws']==b['draws']
            x,y=lineages[0][parent],lineages[1][parent]
            u,v=np.array(x['division_axis']),np.array(y['division_axis'])
            angle=float(np.degrees(np.arccos(np.clip(abs(u@v)/(np.linalg.norm(u)*np.linalg.norm(v)),0,1))))
            records.append(dict(parent=parent,axis_angle_degrees=angle,
                division_start_time_difference=abs(x['division_start']-y['division_start']),
                abscission_time_difference=abs(a['time']-b['time'])))
        rows.append(dict(branch=branch,lineage_and_partition_correspondence=True,
            maximum_axis_angle_degrees=max(r['axis_angle_degrees'] for r in records),parents=records))
    pins[str(Path(__file__).resolve())]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'cleavage-clock-diagnosis.json').write_text(json.dumps({'branches':rows,'source_sha256':pins,
        'interpretation':'Descriptive raw same-parent differences; no causal attribution or time/spatial alignment.'},indent=2)+'\n')


if __name__=='__main__':main()
