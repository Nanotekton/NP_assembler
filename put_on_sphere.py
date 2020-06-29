import numpy as np
import numba 

#@numba.njit
def move_atom(atom, shell, heading_point, sep=3.5):
   '''move point until it touches a convex shell; heading should be the center of the shell'''
   heading = heading_point-atom
   dist = np.linalg.norm(heading)
   assert dist!=0, 'Already there :)'
   heading/=dist
   nearest = None
   nearest_dist = np.inf

   Nshell = shell.shape[0]   

#   for i in numba.prange(Nshell):
   for i in range(Nshell):
      v = shell[i]-atom
      v_m = np.linalg.norm(v)
      vh_m = heading.dot(v)
      vh = heading*vh_m
      vp = v - vh
      vp_m = np.linalg.norm(vp)
      if vp_m>sep: continue
      if v_m<nearest_dist:
         nearest_dist = v_m
         nearest = v
   
   vh_m = heading.dot(nearest)
   vh = heading*vh_m
   vp = nearest - vh
   vp_m = np.linalg.norm(vp)
   neg_sep = (sep**2-vp_m**2)**0.5
   distance_to_travel = vh_m - neg_sep
   
   new_position = atom+heading*distance_to_travel
   return new_position
 
  
#@numba.njit(parallel=True)
def put_all_on_shell(xyz, shell, sep=3.5):
   xyz_center = xyz.mean(axis=0)
   shell_center = shell.mean(axis=0)
   xyz = xyz - xyz_center
   
   rad_x = ((xyz**2).sum(axis=1)**0.5).min()
   rad_s = (((shell-shell_center)**2).sum(axis=1)**0.5).max()
   factor = 1.5*rad_s/rad_x
   
   xyz = xyz*factor + shell_center
   N=xyz.shape[0]
   
#   for i in numba.prange(N):
   for i in range(N):
      try:
         xyz[i] = move_atom(xyz[i], shell, shell_center, sep)
      except:
         print i
         raise
   
   return xyz
   

if __name__=='__main__':
   import argparse
   parser=argparse.ArgumentParser()
   parser.add_argument('--shell', type=str)
   parser.add_argument('--r', type=float, default=1.8)
   parser.add_argument('--r_shell', type=float, default=1.7)
   parser.add_argument('--out', type=str)
   parser.add_argument('input_file', type=str)
   args = parser.parse_args()

   sep = args.r+args.r_shell
   shell = np.loadtxt(args.shell, skiprows=2, usecols=(1,2,3))
   atoms = np.loadtxt(args.input_file, skiprows=2, dtype=str)
   elements = atoms[:,0]
   xyz = atoms[:,1:].astype(float)
   
   new_xyz = put_all_on_shell(xyz, shell, sep)
   with open(args.out, 'w') as f:
      f.write('%5i\n\n'%(len(elements)))
      for i in range(len(elements)):
         f.write('%5s '%elements[i])
         for x in new_xyz[i]:
            f.write('%8.3f '%x)
         f.write('\n')
         

