"""C6 H-matrix: ACA far blocks, streamed full residual certification, direct near blocks.

The FULL residual check makes setup potentially O(Q^2); only repeated
applications can be subquadratic. Local Frobenius tolerances are not a
relative current/field error guarantee. No physics kernel is reimplemented.
"""
from dataclasses import dataclass, asdict
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import numpy as np

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Settings:
    tol: float = 1e-6
    leaf: int = 16
    max_rank: int = 32
    eta: float = 2.0
    phase_limit: float = 8.0
    memory_mib: int = 512

    def validate(self):
        if not (0 < self.tol < 1 and self.leaf >= 2 and self.max_rank >= 1
                and 0 < self.eta <= 100 and 0 < self.phase_limit <= 1e6
                and self.memory_mib >= 16):
            raise ValueError('Invalid H-matrix settings')


class Kernel:
    """Fortran owns one geometry per process; do not interleave configurations."""
    _generation = 0
    def __init__(self, config):
        Kernel._generation += 1
        self.generation = Kernel._generation
        self.lib = C.CDLL(str(ROOT/'build/libhm.so'))
        self.lib.hm_init.argtypes = [C.c_char_p, C.c_int, C.POINTER(C.c_int),
                                      C.POINTER(C.c_int), C.POINTER(C.c_double)]
        self.lib.hm_init.restype = None
        self.lib.hm_coords.argtypes = [C.c_void_p]
        self.lib.hm_coords.restype = None
        self.lib.hm_entries.argtypes = [C.c_int,C.c_int,C.c_void_p,C.c_void_p,C.c_void_p]
        self.lib.hm_entries.restype = None
        self.lib.hm_direct.argtypes = [C.c_int,C.c_int,C.c_int,C.c_void_p,C.c_void_p,C.c_void_p,C.c_void_p]
        self.lib.hm_direct.restype = None
        self.lib.hm_near.argtypes = [C.c_int]*4+[C.c_void_p]*6
        self.lib.hm_near.restype = None
        self.lib.hm_panel_matrix.argtypes = [C.c_int]*3+[C.c_void_p]*3
        self.lib.hm_panel_matrix.restype = None
        self.lib.hm_apply_direct.argtypes = [C.c_void_p,C.c_void_p]
        self.lib.hm_apply_direct.restype = None
        q,qc,k = C.c_int(),C.c_int(),C.c_double()
        path = os.fsencode(Path(config).resolve())
        self.lib.hm_init(path,len(path),C.byref(q),C.byref(qc),C.byref(k))
        self.q,self.qcone,self.k = q.value,qc.value,k.value
        self.xyz = np.empty((6,self.q,3),dtype=np.float64)
        self.lib.hm_coords(self.xyz.ctypes.data)
        self.entries_evaluated = 0

    def entries(self,d,rows,cols):
        if self.generation != Kernel._generation:
            raise RuntimeError('Fortran geometry was replaced by another Kernel instance')
        rows=np.ascontiguousarray(rows,dtype=np.int32)
        cols=np.ascontiguousarray(cols,dtype=np.int32)
        if rows.shape != cols.shape or rows.ndim != 1:
            raise ValueError('Entry index shape mismatch')
        if (np.any(rows<0) or np.any(cols<0) or np.any(rows>=3*self.q)
                or np.any(cols>=3*self.q) or not 0<=d<6):
            raise ValueError('Entry index outside geometry')
        out=np.empty(len(rows),dtype=np.complex128)
        self.lib.hm_entries(d,len(rows),rows.ctypes.data,cols.ctypes.data,out.ctypes.data)
        self.entries_evaluated += len(rows)
        if not np.isfinite(out).all():
            raise RuntimeError('Nonfinite kernel')
        return out

    def apply_direct(self,j):
        if self.generation != Kernel._generation:
            raise RuntimeError('Fortran geometry was replaced by another Kernel instance')
        j=np.ascontiguousarray(j,dtype=np.complex128)
        out=np.empty_like(j)
        self.lib.hm_apply_direct(j.ctypes.data,out.ctypes.data)
        return out

    def matrix(self,d,rows,cols):
        if self.generation != Kernel._generation:
            raise RuntimeError('Fortran geometry was replaced by another Kernel instance')
        if (len(rows)%3==0 and len(cols)%3==0 and
            np.array_equal(rows,dofs(rows[::3]//3)) and
            np.array_equal(cols,dofs(cols[::3]//3))):
            rr=np.ascontiguousarray(rows[::3]//3,dtype=np.int32)
            cc=np.ascontiguousarray(cols[::3]//3,dtype=np.int32)
            out=np.empty((len(rows),len(cols)),complex)
            self.lib.hm_panel_matrix(d,len(rr),len(cc),rr.ctypes.data,cc.ctypes.data,out.ctypes.data)
            self.entries_evaluated+=out.size
            if not np.isfinite(out).all():
                raise RuntimeError('Nonfinite kernel')
            return out
        return self.entries(d,np.repeat(rows,len(cols)),np.tile(cols,len(rows))).reshape(len(rows),len(cols))

    def near(self,d,plan,x):
        if self.generation != Kernel._generation:
            raise RuntimeError('Fortran geometry was replaced by another Kernel instance')
        rowptr,refs,colptr,cols=plan
        x=np.asfortranarray(x,dtype=np.complex128)
        out=np.empty_like(x,order='F')
        self.lib.hm_near(d,len(refs),len(colptr)-1,len(cols),rowptr.ctypes.data,refs.ctypes.data,
                         colptr.ctypes.data,cols.ctypes.data,x.ctypes.data,out.ctypes.data)
        return out

    def direct(self,d,rows,cols,x):
        if self.generation != Kernel._generation:
            raise RuntimeError('Fortran geometry was replaced by another Kernel instance')
        # x is scalar source DOFs x six modes. Fortran arrays are (3,n,6).
        rows=np.ascontiguousarray(rows,dtype=np.int32)
        cols=np.ascontiguousarray(cols,dtype=np.int32)
        x=np.asfortranarray(x,dtype=np.complex128)
        out=np.empty((3*len(rows),6),dtype=np.complex128,order='F')
        self.lib.hm_direct(d,len(rows),len(cols),rows.ctypes.data,cols.ctypes.data,x.ctypes.data,out.ctypes.data)
        return out


def dofs(ids):
    return (3*ids[:,None]+np.arange(3)).ravel().astype(np.int32)


def aca(nr,nc,row,column,tol,max_rank):
    """Partially pivoted complex ACA, A ~= U @ V (no conjugation).

    This only proposes factors. The caller must certify every entry.
    Never interpret a few zero sampled rows as proof of a zero block.
    """
    u=np.empty((nr,0),complex); v=np.empty((0,nc),complex)
    tried=set()
    candidates=np.unique(np.linspace(0,nr-1,min(8,nr),dtype=int))
    i=0
    for rank in range(min(max_rank,nr,nc)):
        residual=row(i)-u[i]@v
        if rank==0 or np.max(abs(residual))<1e-28:
            best=float(np.max(abs(residual)))
            for trial in candidates:
                if int(trial) in tried:
                    continue
                rr=row(int(trial))-u[trial]@v
                if np.max(abs(rr))>best:
                    i=int(trial); residual=rr; best=float(np.max(abs(rr)))
            if best<1e-28:
                return None if rank==0 else (u,v)
        j=int(np.argmax(abs(residual)))
        col=column(j)-u@v[:,j]
        pivot=residual[j]
        if abs(pivot)<1e-28:
            return None
        unew=col/pivot
        u=np.column_stack((u,unew)); v=np.vstack((v,residual))
        tried.add(i)
        # Frobenius norm of U V via small Gram matrices.
        norm2=float(np.real(np.trace((u.conj().T@u)@(v@v.conj().T))))
        increment=float(np.linalg.norm(unew)*np.linalg.norm(residual))
        if rank>=1 and increment<=0.1*tol*np.sqrt(max(norm2,1e-300)):
            return u,v
        magnitudes=abs(unew).copy()
        magnitudes[list(tried)]=0
        i=int(np.argmax(magnitudes))
        if magnitudes[i]==0:
            return u,v
    return u,v


def certify(nr,nc,matrix,u,v,tol,chunk=9):
    """Exhaustive streamed Frobenius residual, not sampled validation."""
    norm2=error2=0.0
    for start in range(0,nr,chunk):
        stop=min(start+chunk,nr)
        a=matrix(start,stop)
        e=a-u[start:stop]@v
        norm2+=float(np.vdot(a,a).real)
        error2+=float(np.vdot(e,e).real)
    valid=np.isfinite(norm2+error2) and error2<=tol*tol*norm2
    return valid,norm2,error2


class Node:
    def __init__(self,ids,kernel,leaf,piece):
        self.ids=np.asarray(ids,dtype=np.int32)
        self.piece=piece
        xyz=kernel.xyz[:,self.ids,:]
        self.lo=xyz.min(axis=1); self.hi=xyz.max(axis=1)
        self.diameter=float(np.linalg.norm(self.hi[0]-self.lo[0]))
        self.children=None
        if len(ids)>leaf:
            axis=int(np.argmax(self.hi[0]-self.lo[0]))
            order=np.argsort(xyz[0,:,axis],kind='stable')
            mid=len(ids)//2
            self.children=(Node(self.ids[order[:mid]],kernel,leaf,piece),
                           Node(self.ids[order[mid:]],kernel,leaf,piece))


class HMatrix:
    def __init__(self,kernel,settings,cache,identity):
        settings.validate()
        self.kernel,self.settings = kernel,settings
        self.cache=Path(cache); self.identity=identity
        self.blocks=[]; self.bytes=0
        self.stats=dict(low_rank_blocks=0,direct_blocks=0,zero_blocks=0,
                        certified_norm2=0.0,certified_error2=0.0,stored_bytes=0,
                        cache_loaded=False)
        self._last_progress=time.monotonic()
        if (self.cache/'manifest.json').exists():
            self._load()
        else:
            t=time.monotonic()
            roots=[Node(np.arange(0,kernel.qcone),kernel,settings.leaf,0),
                   Node(np.arange(kernel.qcone,kernel.q),kernel,settings.leaf,1)]
            for d in range(6):
                for a in roots:
                    for b in roots:
                        self._build(d,a,b)
            self.stats.update(build_seconds=time.monotonic()-t,
                              scalar_kernel_evaluations=kernel.entries_evaluated)
            self._save()
        self._prepare_near()
        print('H-matrix '+json.dumps(self.stats),flush=True)

    def _append(self,d,a,b,u=None,v=None,norm2=0,error2=0):
        rows=a.ids.copy(); cols=b.ids.copy()
        needed=rows.nbytes+cols.nbytes+(0 if u is None else u.nbytes+v.nbytes)
        # Include conservative Python/container overhead per block.
        self.bytes+=needed+1024
        if self.bytes>self.settings.memory_mib*2**20:
            raise MemoryError('H-matrix cache budget exceeded; previous order retained')
        if len(self.blocks)>=200000:
            raise MemoryError('H-matrix block limit exceeded; previous order retained')
        self.blocks.append((d,rows,cols,u,v))
        self.stats['direct_blocks' if u is None else 'low_rank_blocks']+=1
        self.stats['certified_norm2']+=norm2
        self.stats['certified_error2']+=error2
        self.stats['stored_bytes']=self.bytes

    def _build(self,d,a,b):
        if d==0 and a.piece==b.piece:
            self.stats['zero_blocks']+=1
            return
        gap=np.maximum(0,np.maximum(a.lo[0]-b.hi[d],b.lo[d]-a.hi[0]))
        distance=float(np.linalg.norm(gap))
        admissible=(distance>0 and max(a.diameter,b.diameter)<=self.settings.eta*distance
                    and self.kernel.k*a.diameter*b.diameter<=self.settings.phase_limit*distance)
        nr,nc=3*len(a.ids),3*len(b.ids)
        scratch=16*(nr+nc)*(self.settings.max_rank+16)
        if admissible and min(nr,nc)>6 and scratch<self.settings.memory_mib*2**19:
            rows,cols=dofs(a.ids),dofs(b.ids)
            factors=aca(nr,nc,
                lambda i:self.kernel.matrix(d,rows[i:i+1],cols)[0],
                lambda j:self.kernel.matrix(d,rows,cols[j:j+1])[:,0],
                self.settings.tol,self.settings.max_rank)
            if factors is not None:
                u,v=factors
                if u.size+v.size<0.8*nr*nc:
                    valid,norm2,error2=certify(nr,nc,
                        lambda start,stop:self.kernel.matrix(d,rows[start:stop],cols),
                        u,v,self.settings.tol)
                    if valid:
                        self._append(d,a,b,u,v,norm2,error2)
                        return
        if a.children is None and b.children is None:
            self._append(d,a,b)
        elif b.children is None or (a.children is not None and a.diameter>=b.diameter):
            for child in a.children:
                self._build(d,child,b)
        else:
            for child in b.children:
                self._build(d,a,child)
        if time.monotonic()-self._last_progress>10:
            print(f'H-matrix building: {len(self.blocks)} blocks, {self.bytes/2**20:.1f} MiB',flush=True)
            self._last_progress=time.monotonic()

    def _prepare_near(self):
        self.near_plans={}
        if self.stats['low_rank_blocks']==0:
            return
        extra=0
        for d in range(6):
            blocks=[(r,c) for dd,r,c,u,v in self.blocks if dd==d and u is None]
            if not blocks: continue
            refs=[[] for _ in range(self.kernel.q)]
            colptr=[0]; columns=[]
            for i,(rows,cols) in enumerate(blocks):
                columns.append(cols)
                colptr.append(colptr[-1]+len(cols))
                for row in rows:
                    refs[row].append(i)
            rowptr=np.concatenate(([0],np.cumsum([len(x) for x in refs])))
            if max(rowptr[-1],colptr[-1])>=2**31:
                raise MemoryError('Near plan exceeds 32-bit indexing')
            plan=(rowptr.astype(np.int32),np.array([b for rr in refs for b in rr],dtype=np.int32),
                  np.array(colptr,dtype=np.int32),np.concatenate(columns))
            extra+=sum(x.nbytes for x in plan)
            if self.bytes+extra>self.settings.memory_mib*2**20:
                raise MemoryError('Near plan exceeds cache memory budget')
            self.near_plans[d]=plan
        self.stats['near_plan_bytes']=extra

    def apply(self,j):
        q=self.kernel.q
        if np.asarray(j).size!=18*q or not np.isfinite(j).all():
            raise ValueError('Invalid current')
        if self.stats['low_rank_blocks']==0:
            return self.kernel.apply_direct(np.asarray(j).ravel())
        modes=np.fft.fft(np.asarray(j).reshape(6,q,3),axis=0)
        x=modes.reshape(6,3*q).T
        y=np.zeros_like(x)
        phases=np.exp(2j*np.pi*np.arange(6)[:,None]*np.arange(6)/6)
        for d,plan in self.near_plans.items():
            y+=self.kernel.near(d,plan,x)*phases[d]
        for d,rows,cols,u,v in self.blocks:
            if u is None: continue
            rd,cd=dofs(rows),dofs(cols)
            source=x[cd]
            block=u@(v@source)
            y[rd]+=block*phases[d]
        result=np.fft.ifft(y.T.reshape(6,q,3),axis=0).ravel()
        if not np.isfinite(result).all():
            raise RuntimeError('Nonfinite H-matrix current; previous checkpoint retained')
        return result

    def _save(self):
        # The run output lock also protects its private cache.
        self.cache.mkdir(parents=True,exist_ok=True)
        arrays={}
        for i,(d,r,c,u,v) in enumerate(self.blocks):
            arrays[f'd{i}']=np.array(d,dtype=np.int32)
            arrays[f'r{i}']=r; arrays[f'c{i}']=c
            if u is not None:
                arrays[f'u{i}']=u; arrays[f'v{i}']=v
        fd,tmp=tempfile.mkstemp(prefix='.pending-',dir=self.cache)
        with os.fdopen(fd,'wb') as f:
            np.savez(f,**arrays); f.flush(); os.fsync(f.fileno())
        checksum=sha256(Path(tmp))
        os.replace(tmp,self.cache/'blocks.npz')
        manifest=dict(schema=1,identity=self.identity,settings=asdict(self.settings),
                      count=len(self.blocks),sha256=checksum,stats=self.stats)
        fd,tmp=tempfile.mkstemp(prefix='.pending-',dir=self.cache)
        with os.fdopen(fd,'w') as f:
            json.dump(manifest,f,allow_nan=False); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,self.cache/'manifest.json')
        fd=os.open(self.cache,os.O_RDONLY|os.O_DIRECTORY)
        try: os.fsync(fd)
        finally: os.close(fd)

    def _load(self):
        meta=json.loads((self.cache/'manifest.json').read_text())
        if (meta['schema']!=1 or meta['identity']!=self.identity
                or meta['settings']!=asdict(self.settings)):
            raise RuntimeError('H-matrix cache identity mismatch')
        path=self.cache/'blocks.npz'
        if path.stat().st_size>self.settings.memory_mib*2**20*2:
            raise MemoryError('Cache file exceeds budget')
        if sha256(path)!=meta['sha256']:
            raise RuntimeError('H-matrix cache checksum mismatch; cache preserved')
        with np.load(path,allow_pickle=False) as data:
            expected=set()
            for i in range(meta['count']):
                expected.update([f'd{i}',f'r{i}',f'c{i}'])
                if f'u{i}' in data:
                    expected.update([f'u{i}',f'v{i}'])
            if set(data.files)!=expected:
                raise RuntimeError('H-matrix cache manifest/array mismatch')
            for i in range(meta['count']):
                d=int(data[f'd{i}']); r=data[f'r{i}']; c=data[f'c{i}']
                u=data[f'u{i}'] if f'u{i}' in data else None
                v=data[f'v{i}'] if u is not None else None
                if (not 0<=d<6 or r.ndim!=1 or c.ndim!=1 or
                    r.dtype!=np.int32 or c.dtype!=np.int32 or
                    len(r)==0 or len(c)==0 or min(r.min(),c.min())<0 or
                    max(r.max(),c.max())>=self.kernel.q):
                    raise RuntimeError('Invalid H-matrix cache indices')
                if u is not None and (u.ndim!=2 or v.ndim!=2 or
                    u.shape[0]!=3*len(r) or v.shape!=(u.shape[1],3*len(c)) or
                    not np.isfinite(u).all() or not np.isfinite(v).all()):
                    raise RuntimeError('Invalid H-matrix cache factors')
                self.blocks.append((d,r,c,u,v))
                self.bytes+=r.nbytes+c.nbytes+1024+(0 if u is None else u.nbytes+v.nbytes)
                if self.bytes>self.settings.memory_mib*2**20:
                    raise MemoryError('Loaded cache exceeds budget')
        self.stats=meta['stats']
        low=sum(u is not None for _,_,_,u,_ in self.blocks)
        if low!=self.stats['low_rank_blocks'] or len(self.blocks)-low!=self.stats['direct_blocks']:
            raise RuntimeError('H-matrix cache block count mismatch')
        self.stats['cache_loaded']=True


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()
