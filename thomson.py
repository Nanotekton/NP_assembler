import numpy as np
from sys import stdout
import numba

@numba.njit(parallel=True)
def compute_sep_tensor(xyz):
   N=xyz.shape[0]
   dist = np.zeros((N, N, 3))
   for i in numba.prange(N):
      for j in numba.prange(i+1, N):
         d=xyz[j]-xyz[i]
         dist[i,j] = d
         dist[j,i] = -d
   return dist


def compute_charge_mask(qs):
   qs = qs.reshape(-1,1)
   return qs.dot(qs.T)


@numba.njit(parallel=True)
def calc_force(dist_tensor, charge_mask):
   dists = (dist_tensor**2).sum(axis=2)
   N = dists.shape[0]
   factor = np.zeros((N,N))
   for i in numba.prange(N):
      for j in numba.prange(i+1, N):
         factor[i,j] = charge_mask[i,j]/dists[i,j]**1.5
         factor[j,i] = factor[i,j]
   factor = factor.reshape(N,N,1)
   f = (factor*dist_tensor).sum(axis=0)
   return f


def normalize(xyz):
   mod = np.sqrt((xyz*xyz).sum(axis=1).reshape(-1,1))
   factor = np.where(mod!=0,1/mod,0)
   return xyz*factor, mod


def remove_radial_component(xyz, force):
   '''CENTER ASSUMED TO BE IN 0,0,0 '''
   rs, _ = normalize(xyz)
   scalar_prod = (force*rs).sum(axis=1)
   radial = scalar_prod.reshape(-1,1)*rs
   return force-radial


def pull_to_the_sphere(xyz, radius):
   normalized, _ = normalize(xyz)
   return normalized*radius 


def thomson(qs, crd, tol=4.5e-4, maxit=10000, stepsize=0.01, radius=22, v=True, center=np.array([0,0,0]), fixed=None):
   xyz = crd-center
   charge_mask = compute_charge_mask(qs)
   max_f = 1000
   vectors = compute_sep_tensor(xyz)
   iteration=0
   while iteration<maxit and max_f>tol:
      force = calc_force(vectors, charge_mask)
      force = remove_radial_component(xyz, force)
      force, magnitudes = normalize(force)
      if type(fixed).__name__!='NoneType':
         force = force*(1-fixed.reshape(-1,1))
         magnitudes = magnitudes*(1-fixed)
      magnitudes = np.where(magnitudes>1,1,magnitudes)
      max_f = magnitudes.max()
      xyz += stepsize*force
      xyz = pull_to_the_sphere(xyz, radius)
      vectors = compute_sep_tensor(xyz)
      iteration+=1
      if v:
         stdout.write('Iter: %i  Max F: %10.6f\r'%(iteration, max_f))
         stdout.flush()
   if iteration>maxit:print('\nmaxit exceeded')
   if v:
         stdout.write('\nDone\n')
         stdout.flush()
   return xyz+center


def gen_attachement_points(qs, xyz, enhance_half=False, fix_half=False, pull=False, maxit=1000, radius=22, factor=0.5):
   qs=-np.ones(qs.shape)
   if enhance_half:
      indices = np.arange(qs.shape[0])
      half = np.random.choice(indices, int(indices.shape[0]*factor), False)
      qs[half]=-2
   else:
      qs=-np.ones(qs.shape)
   if pull:
      xyz-=xyz.mean(axis=0)
      xyz=pull_to_the_sphere(xyz, radius)
   if fix_half:
      if enhance_half:
         which_first =np.where(qs==-1)[0]
      else:
         which_first = np.arange(int(len(qs)/2))
      qs_first, xyz_first = qs[which_first], xyz[which_first,:]
      xyz[which_first,:] = thomson(qs_first, xyz_first, maxit=maxit)
      fix_mask = np.zeros(len(qs))
      fix_mask[which_first] = 1.0
   else:
      fix_mask = None
   new_xyz = thomson(qs, xyz, maxit=maxit, fixed=fix_mask)
   
   return qs, new_xyz



if __name__=='__main__':
   import argparse
   parser=argparse.ArgumentParser()
   parser.add_argument('--out', default='thomson_opt.xyz', type=str)
   parser.add_argument('--pull', action='store_true')
   parser.add_argument('--enhance_half', action='store_true')
   parser.add_argument('--fix_half', action='store_true')
   parser.add_argument('--maxit', type=int, default=10000)
   parser.add_argument('--center', type=float, nargs=3, default=[0,0,0])
   parser.add_argument('input_file', type=str)
   args = parser.parse_args()
   
   data = np.loadtxt(args.input_file)
   qs, xyz = data[:,0], data[:,1:]
   qs, new_xyz  = gen_attachement_points(qs, xyz, enhance_half=args.enhance_half, fix_half=args.fix_half, pull=args.pull, maxit=args.maxit, radius=22, factor=0.5)
   with open(args.out, 'w') as f:
      f.write('%i\n\n'%len(qs))
      for i,x in enumerate(qs):
         if x==-1: s='O'
         elif x==-2: s='N'
         else: s='H'
         f.write('%5s %8.3f %8.3f %8.3f\n'%(s, new_xyz[i,0], new_xyz[i,1], new_xyz[i,2]))
   print('\nDone')
   
