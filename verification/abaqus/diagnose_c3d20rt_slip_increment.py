"""Compare finite slip increment definitions at prescribed native geometries.

Uses the nearest primary projection, with edge extrapolation, as an approximation; it
does not solve equilibrium or substitute for the production disk transfer.
"""
import numpy as np
import netCDF4 as nc
import diagnose_c3d20rt_finite_friction as d

def shapes(x, y):
    sx = np.array([-1, 1, 1, -1]); sy = np.array([-1, -1, 1, 1])
    x = x[:, None]; y = y[:, None]
    n = np.concatenate([.25*(1+sx*x)*(1+sy*y)*(sx*x+sy*y-1),
        .5*(1-x*x)*(1-y), .5*(1+x)*(1-y*y), .5*(1-x*x)*(1+y), .5*(1-x)*(1-y*y)], axis=1)
    dx = np.concatenate([.25*sx*(1+sy*y)*(2*sx*x+sy*y),
        -x*(1-y), .5*(1-y*y), -x*(1+y), -.5*(1-y*y)], axis=1)
    dy = np.concatenate([.25*sy*(1+sx*x)*(sx*x+2*sy*y),
        -.5*(1-x*x), -(1+x)*y, .5*(1-x*x), -(1-x)*y], axis=1)
    return n, dx, dy

def project(p, face):
    xy = np.zeros((len(p), 2))
    for _ in range(80):
        n, dx, dy = shapes(*xy.T)
        q = np.einsum('ti,tic->tc', n, face)
        A = np.einsum('kti,tic->tkc', np.array([dx, dy]), face)
        rhs = np.einsum('tkc,tc->tk', A, p-q)
        delta = np.linalg.solve(A@A.transpose(0,2,1), rhs[...,None])[...,0]
        xy += delta
        valid = np.max(abs(xy),axis=1) <= 1+1e-8
        if not valid.any() or abs(delta[valid]).max() < 1e-13: break
    else: raise RuntimeError('Valid projection did not converge')
    n = shapes(*xy.T)[0]
    q = np.einsum('ti,tic->tc', n, face)
    dist = np.linalg.norm(p-q, axis=1)
    dist += np.linalg.norm(np.maximum(abs(xy)-1,0),axis=1)
    return n, q, dist

def basis(face, dx, dy):
    tx = np.einsum('i,tic->tc', dx, face); ty = np.einsum('i,tic->tc', dy, face)
    av = -np.cross(tx, ty); measure = np.linalg.norm(av, axis=1)
    normal = av/measure[:,None]
    first = ty/np.linalg.norm(ty,axis=1)[:,None]
    return np.array([first,np.cross(normal,first)]).transpose(1,0,2), normal, measure, np.array([tx,ty]).transpose(1,0,2)

with nc.Dataset(d.ROOT/'verification/meshes/c3d20rt_finite_contact.e') as ds:
    i = nc.chartostring(ds['ss_names'][:]).tolist().index('primary_contact')+1
    primaryids = [d.elements[int(e)-1,d.FACES[int(s)-1]] for e,s in zip(ds['elem_ss'+str(i)][:],ds['side_ss'+str(i)][:])]
old = np.concatenate([d.x0[None], d.X[:-1]])
increment = {key:np.zeros((len(d.X),13,2)) for key in ['current_virtual','local_secant','current_projection','midpoint_projection','constraint_secant','old_primary_virtual','mid_primary_virtual','mean_basis_virtual','corotated_virtual']}
reference_increment = np.zeros_like(increment['current_virtual'])
reference_area = np.zeros((len(d.X),13))
mid_increment = np.zeros_like(reference_increment)
mid_area = np.zeros_like(reference_area)
area = np.zeros((len(d.X),13))
old_mean = np.concatenate([np.broadcast_to(np.array([[0,0,1],[0,-1,0]]),(1,13,2,3)), d.basis[:-1]])
for ids in d.faces:
    face = d.X[:,ids]; old_face = old[:,ids]
    for local, node in enumerate(ids):
        row = d.lookup[int(node)]
        for a,b,w in d.samples(local):
            n,dx,dy = d.shape(2*a-1,2*b-1)
            p = np.einsum('i,tic->tc',n,face)
            po = np.einsum('i,tic->tc',n,old_face)
            current_basis,normal,measure,A = basis(face,dx,dy)
            old_basis,_,_,_ = basis(old_face,dx,dy)
            mid_basis,_,mid_measure,_ = basis((face+old_face)/2,dx,dy)
            _,_,ref_measure,_ = basis(d.x0[None,ids],dx,dy)
            projections = [project(p,d.X[:,pi]) for pi in primaryids]
            old_projections = [project(po,old[:,pi]) for pi in primaryids]
            old_selected = np.argmin(np.array([v[2] for v in old_projections]),axis=0)
            selected = np.argmin(np.array([v[2] for v in projections]),axis=0)
            if not np.isfinite(np.min(np.array([v[2] for v in projections]),axis=0)).all():
                raise RuntimeError('No valid primary projection')
            q = np.zeros_like(p); qo = np.zeros_like(p)
            for pi,(pn,pq,dist),index in zip(primaryids,projections,range(len(primaryids))):
                use = selected==index
                q[use] = pq[use]
                qo[use] = np.einsum('ti,tic->tc',pn,old[:,pi])[use]
            relative = q-p; old_relative = qo-po
            du = old_relative-relative
            old_primary_du = np.zeros_like(p)
            for pi,(pn,_,_),index in zip(primaryids,old_projections,range(len(primaryids))):
                use = old_selected==index
                old_primary_du[use] = np.einsum('ti,tic->tc',pn,d.X[:,pi]-old[:,pi])[use]
            old_du = p-po-old_primary_du
            gap = -np.sum(relative*normal,axis=1)
            gradient = np.array([dx,dy]).T[None]@np.linalg.solve(A@A.transpose(0,2,1),A)
            rotation = np.einsum('tic,tkc,tij,tj->tk',gradient,current_basis,face-old_face,normal)
            wt = 4*(1 if local<4 else 5)/24*w*measure
            values = {
                'current_projection':np.einsum('tc,tkc->tk',du,current_basis),
                'midpoint_projection':np.einsum('tc,tkc->tk',du,mid_basis),
                'local_secant':np.einsum('tc,tkc->tk',old_relative,old_basis)-np.einsum('tc,tkc->tk',relative,current_basis),
                'constraint_secant':np.einsum('tc,tkc->tk',old_relative,old_mean[:,row])-np.einsum('tc,tkc->tk',relative,d.basis[:,row]),
            }
            values['current_virtual'] = values['current_projection']+gap[:,None]*rotation
            values['old_primary_virtual'] = np.einsum('tc,tkc->tk',old_du,current_basis)+gap[:,None]*rotation
            values['mid_primary_virtual'] = (values['old_primary_virtual']+values['current_virtual'])/2
            mean_rotation = np.einsum('tic,tkc,tij,tj->tk',gradient,d.basis[:,row],face-old_face,normal)
            values['mean_basis_virtual'] = np.einsum('tc,tkc->tk',du,d.basis[:,row])+gap[:,None]*mean_rotation
            cb = np.concatenate([d.basis[:,row],np.cross(d.basis[:,row,0],d.basis[:,row,1])[:,None]],axis=1)
            ob = np.concatenate([old_mean[:,row],np.cross(old_mean[:,row,0],old_mean[:,row,1])[:,None]],axis=1)
            R = cb.transpose(0,2,1)@ob
            corotated_du = np.einsum('tij,tj->ti',R,old_relative)-relative
            corotated_face_du = face-np.einsum('tij,tkj->tki',R,old_face)
            corotated_rotation = np.einsum('tic,tkc,tij,tj->tk',gradient,current_basis,corotated_face_du,normal)
            values['corotated_virtual'] = np.einsum('tc,tkc->tk',corotated_du,current_basis)+gap[:,None]*corotated_rotation
            wr = 4*(1 if local<4 else 5)/24*w*ref_measure
            wm = 4*(1 if local<4 else 5)/24*w*mid_measure
            reference_increment[:,row] += wr[:,None]*values['current_virtual']
            reference_area[:,row] += wr
            mid_increment[:,row] += wm[:,None]*values['current_virtual']
            mid_area[:,row] += wm
            for key in increment: increment[key][:,row] += wt[:,None]*values[key]
            area[:,row] += wt
native = d.frames[:,:,8:10]
native_increment = np.diff(np.concatenate([np.zeros_like(native[:1]),native]),axis=0)
increment['reference_weight_virtual'] = reference_increment/reference_area[:,:,None]*area[:,:,None]
increment['midpoint_weight_virtual'] = mid_increment/mid_area[:,:,None]*area[:,:,None]
for key,value in increment.items():
    value /= area[:,:,None]
    difference = np.cumsum(value,axis=0)-native
    print(key,'increment_maximum_absolute_m',np.linalg.norm(value-native_increment,axis=2).max(),
          'accumulated_relative_l2',np.linalg.norm(difference)/np.linalg.norm(native),
          'accumulated_maximum_absolute_m',np.linalg.norm(difference,axis=2).max(),flush=True)

