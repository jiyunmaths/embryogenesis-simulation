"""Visualization-only isosurfaces of the actual phase field; no smoothing of phi."""
from pathlib import Path

import numpy as np
from skimage.measure import marching_cubes


def cell_mesh(field, dx, extent, level=.5):
    """Extract a cropped, full-resolution triangle mesh in world coordinates.

    No padding, convex hull, artificial cap, or shape fitting is used. Surfaces
    intersecting the computational boundary remain open and are labeled as such.
    """
    field=np.asarray(field)
    if field.ndim!=3 or not np.isfinite(field).all():
        raise ValueError('mesh input must be a finite 3D field')
    empty={'vertices':[],'faces':[],'normals':[],'closed':False,'boundary_edges':0,'level':level}
    inside=np.argwhere(field>level)
    if not len(inside):return empty
    lo=np.maximum(inside.min(axis=0)-1,0)
    hi=np.minimum(inside.max(axis=0)+2,field.shape)
    crop=field[tuple(slice(a,b) for a,b in zip(lo,hi))]
    if min(crop.shape)<2 or crop.min()>=level:return empty
    vertices,faces,normals,_=marching_cubes(crop,level=level,spacing=(dx,dx,dx),
                                          method='lewiner',allow_degenerate=False)
    vertices+=((lo+.5)*dx-extent)[None,:]
    # Align triangle winding with the outward gradient normals, allowing correct
    # back-face lighting without changing the extracted surface geometry.
    triangles=vertices[faces]
    face_normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    if np.sum(face_normals*normals[faces].mean(axis=1))<0:
        faces=faces[:,::-1].copy()
        face_normals=-face_normals
    # Smooth lighting normals on the extracted geometry. No vertex is moved;
    # averaging adjacent face areas avoids grid-gradient lighting artifacts.
    normals=np.zeros_like(vertices)
    for corner in range(3):
        np.add.at(normals,faces[:,corner],face_normals)
    normals/=np.maximum(np.linalg.norm(normals,axis=1,keepdims=True),1e-12)
    edges=np.sort(np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]])),axis=1)
    _,counts=np.unique(edges,axis=0,return_counts=True)
    boundary=int(np.sum(counts==1))
    return {'vertices':vertices.tolist(),'faces':faces.tolist(),'normals':normals.tolist(),
            'closed':bool(len(faces) and np.all(counts==2)), 'boundary_edges':boundary,'level':level}


def viewer_template():
    """Inline the shared renderer so saved viewers work without a server/CDN."""
    directory=Path(__file__).parent
    return (directory/'viewer.html').read_text().replace(
        '__CELL_SURFACE_RENDERER__',(directory/'cell_surface.js').read_text())


def main():
    """Export an exact checkpoint snapshot; old point clouds are not remeshed."""
    import argparse
    from dataclasses import asdict
    import hashlib
    import json
    from .model import Simulation
    parser=argparse.ArgumentParser(description='Export shaded cell surfaces from a full simulation checkpoint.')
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('choose a new output directory')
    sim=Simulation.restore(args.checkpoint)
    frame={'metrics':sim.metrics(),'cells':sim.surfaces(),'graph':sim.graph_snapshot()}
    payload=json.dumps({'config':asdict(sim.config),'frames':[frame]},separators=(',',':'),allow_nan=False)
    template=viewer_template().replace('One cell. Two possible identities.','Cell surfaces in 3D.')
    template=template.replace('A freely evolving 3D aggregate of deformable cells.',
                              f'Saved state at simulation time {sim.time:g}. Drag to explore the cell boundaries.')
    args.output.mkdir(parents=True)
    (args.output/'trajectory.json').write_text(payload)
    (args.output/'viewer.html').write_text(template.replace('__SIMULATION_DATA__',payload))
    (args.output/'provenance.json').write_text(json.dumps({'checkpoint':str(args.checkpoint),
        'sha256':hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),'time':sim.time,
        'scope':'One exact saved state; no reconstruction of earlier point-only frames and no solver steps.'},indent=2)+'\n')
    print(f'Wrote {args.output / "viewer.html"}')


if __name__=='__main__':main()
