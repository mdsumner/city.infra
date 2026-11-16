"""Network health checks on the tas-transport-web segments layer.

Mirrors what sfnetworks construction and cleaning would do (node identity,
smoothing, subdivision, snapping, noding at crossings), using the node0/node1
ids already in the layer. Python only because R was not available in the
sandbox; the R version would use wkpool (see sfnetworks-wkpool.md).

Usage:
  curl -o segs.parquet \
    https://data.source.coop/mdsumner/tas-transport-web/2025-07-10/list_transport_segments.parquet
  python network_checks.py segs.parquet

Needs: numpy, scipy, pyarrow, shapely (>= 2).
"""
import sys, numpy as np, pyarrow.parquet as pq, shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
t = pq.read_table(sys.argv[1])
n0 = t['node0'].to_numpy(); n1 = t['node1'].to_numpy()
ne = len(n0)
ids, inv = np.unique(np.concatenate([n0,n1]), return_inverse=True)
i0, i1 = inv[:ne], inv[ne:]
nn = len(ids)
print("edges", ne, "nodes", nn, "id range", ids.min(), ids.max())
loop = i0 == i1
print("self loops", loop.sum())
lo, hi = np.minimum(i0,i1), np.maximum(i0,i1)
pair = lo.astype(np.int64)*nn + hi
u, c = np.unique(pair[~loop], return_counts=True)
print("multi-edge pairs", (c>1).sum(), "extra edges", (c-1).sum())
deg = np.bincount(np.concatenate([i0,i1]), minlength=nn)
print("degree dist", {k:int(v) for k,v in zip(*np.unique(np.minimum(deg,6), return_counts=True))})
A = coo_matrix((np.ones(ne), (i0,i1)), shape=(nn,nn))
ncomp, lab = connected_components(A, directed=False)
sz = np.bincount(lab)
print("components", ncomp, "largest", sz.max(), f"({sz.max()/nn:.3f} of nodes)", "singleton-edge comps", (sz==2).sum())
pseudo = deg == 2
# pure cycle components: all nodes deg 2
allpseudo = np.bincount(lab, weights=~pseudo) == 0
print("pseudo (deg2) nodes", pseudo.sum(), "pure-cycle comps", allpseudo.sum())
print("edges after smoothing (approx)", ne - pseudo.sum() + allpseudo.sum())
# require_equal on TRANS_TYPE / PRI_NAME: pseudo node kept only if incident edges agree
for col in ["TRANS_TYPE","PRI_NAME"]:
    v = np.array(t[col].fill_null("").to_pylist(), dtype=object)
    _, code = np.unique(v, return_inverse=True)
    # for each node collect min and max code of incident edges
    nodes = np.concatenate([i0,i1]); codes = np.concatenate([code,code])
    mn = np.full(nn, 10**9); mx = np.full(nn, -1)
    np.minimum.at(mn, nodes, codes); np.maximum.at(mx, nodes, codes)
    ok = pseudo & (mn == mx)
    print(f"pseudo nodes with equal {col}", ok.sum(), "-> edges after", ne-ok.sum())
# geometry checks
g = shapely.from_wkb(t['geometry'].to_numpy(zero_copy_only=False))
print("geom types", np.unique(shapely.get_type_id(g), return_counts=True))
coords, idx = shapely.get_coordinates(g, return_index=True)
first = np.r_[True, idx[1:] != idx[:-1]]
last = np.r_[idx[1:] != idx[:-1], True]
s = coords[first]; e = coords[last]
# sfnetworks-style node identity: paste(round(coord, 12))
def key(a, d): return np.char.add(np.char.add(np.round(a[:,0],d).astype(str), " "), np.round(a[:,1],d).astype(str))
for d in [12, 7, 6]:
    k = key(np.vstack([s,e]), d)
    _, kid = np.unique(k, return_inverse=True)
    print(f"string-key nodes at {d} digits:", kid.max()+1)
# exact bit-pattern identity
ep = np.vstack([s,e]); _, eid = np.unique(ep.view([('x',float),('y',float)]), return_inverse=True)
eid = eid.ravel()
print("exact-coordinate nodes", eid.max()+1)
# does node id <-> exact coord agree 1:1?
pairs = np.unique(np.stack([inv, eid],1), axis=0)
print("distinct (nodeid, exactcoord) pairs", len(pairs))
# metric coords approx (local equirectangular)
lat0 = np.deg2rad(-42.0)
def m(a): return np.c_[a[:,0]*111320*np.cos(lat0), a[:,1]*110574]
# near-miss endpoints: distinct node ids within 1 m / 0.1 m (undershoot / snapping)
nodexy = np.zeros((nn,2)); nodexy[inv] = m(ep)
tree = cKDTree(nodexy)
for r in [0.01, 0.1, 1.0]:
    pr = tree.query_pairs(r, output_type='ndarray')
    print(f"node pairs closer than {r} m but distinct ids:", len(pr))
# interior vertices coincident with another segment endpoint (vertex T-junction, sfnetworks subdivision case)
interior = ~(first|last)
ic = coords[interior]
icm = m(ic)
d, j = tree.query(icm, distance_upper_bound=1e-6)
hit = np.isfinite(d)
print("interior vertices exactly on a node:", hit.sum(), "segments affected", len(np.unique(idx[interior][hit])))
# dangling deg1 nodes near another segment (undershoot candidates) within 2 m of any vertex not their own
dang = np.where(deg==1)[0]
print("dangling nodes", len(dang))
allm = m(coords); vt = cKDTree(allm)
# which segment does each node belong to (for deg1, exactly one)
seg_of_node = np.full(nn, -1); seg_of_node[i0] = np.arange(ne); seg_of_node[i1] = np.arange(ne)
res = vt.query_ball_point(nodexy[dang], 2.0)
cnt = sum(1 for k,lst in zip(dang,res) if any(idx[q] != seg_of_node[k] for q in lst))
print("dangling nodes within 2 m of another segment's vertex:", cnt)
# segment pairs whose interiors cross with no shared vertex (noding candidates,
# many will be grade separations: bridges, overpasses, rail)
tree2 = shapely.STRtree(g)
a, b = tree2.query(g, predicate="crosses")
k = a < b
print("crossing pairs without a shared vertex:", k.sum())
