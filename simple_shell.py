import argparse
import os, sys
mypath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(mypath)
from v2w import *

line_fmt='HETATM%5i %-4s %3s %5i    %8.3f%8.3f%8.3f  1.00  0.00       %5s\n'

def min_dist(ar1, ar2):
   dmin=100
   for x in ar1:
      for y in ar2:
         d=y-x
         d=d.dot(d)**0.5
         if d<dmin:
            dmin=d
   return dmin

def make_random_z_rotation():
   angle = np.random.random()*np.pi*2
   cos_angle = np.cos(angle)
   sin_angle = np.sin(angle)
   rot = np.array([[cos_angle, -sin_angle, 0], [sin_angle, cos_angle,0], [0,0,1]])
   return rot


def align_along_longest(xyz, center_idx):
   #find largest distance
   N=xyz.shape[0]
   max_dist=0
   dist_idx=-1
   xyz-=xyz[center_idx]
   for i in range(N):
      d=xyz[i]-xyz[center_idx]
      d=d.dot(d)
      if d>max_dist:
         max_dist=d
         dist_idx=i
   assert dist_idx>-1, 'Nothing found'
   r = xyz[dist_idx]
   r = r.dot(r)**0.5+10
   rot, trans = get_rot(xyz[dist_idx], np.array([0,0,r]))
   new = rot.dot(xyz.T).T+trans

   #add random rotation along z-axis
   rot = make_random_z_rotation()
   #new-=trans 
   new = rot.dot(new.T).T#+trans
   
   return new 


def loadXYZ(name):
   data=np.loadtxt(name,skiprows=2,dtype=str)
   names=data[:,0]
   xyz=data[:,1:].astype(float)
   return names,xyz


def write(f,names, names_f,new,res, at, resname):
   for i in range(len(names)):
      at+=1
      at_symb = names_f[i]
      f.write(line_fmt%(at, at_symb, resname, res, new[i][0], new[i][1], new[i][2], names[i]))
   return at


def gen_lines(names, names_f, new, res, at, resname):
   result = []
   for i in range(len(names)):
      at+=1
      at_symb = names_f[i]
      result.append(line_fmt%(at, at_symb, resname, res, new[i][0], new[i][1], new[i][2], names[i]))
   return result, at


def stats_of_two_arrs(r1,r2=None):
   dmin=100
   dav=0
   nr=0
   for i,x in enumerate(r1):
      if r2 is None:
         second = r1[i+1:]
      else:
         second=r2
      for y in second:
         d=min_dist(x,y)
         dav+=d
         nr+=1
         if d< dmin: dmin=d
   dav/=nr
   return dav, dmin


def attach_monomers(monomer_xyz, monomer_anchor_idx, attachement_shell, offset=0):
   monomer_aligned = align_along_longest(monomer_xyz, monomer_anchor_idx)
   anchor = monomer_aligned[monomer_anchor_idx]
   residues = []
   shell_center = attachement_shell.mean(axis=0)
   for point in attachement_shell:
      rot = make_random_z_rotation()
      monomer_aligned = rot.dot(monomer_aligned.T).T
      rot, trans = get_rot(anchor, point-shell_center, offset)
      new = rot.dot(monomer_aligned.T).T+trans+shell_center
      residues.append(new)
   return residues



if __name__=='__main__':    
   parser=argparse.ArgumentParser(description='Make SAM', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--shell',type=str,help='XYZ of outer atoms')
   parser.add_argument('--monomer',type=str,help='XYZ of monomer')
   parser.add_argument('--names',type=str,help='XYZ of monomer')
   parser.add_argument('--n0', type=int, help="index of monomer's atom to be attached", default=0)
   parser.add_argument('--scale', type=float, help="scale factor", default=1.0)
   #parser.add_argument('--d0', type=float, help="minimum distance between monomers", default=2.0)
   parser.add_argument('--out', type=str, help="output file", default='out.pdb')
   
   args=parser.parse_args()
   
   
   with open(args.names, 'r') as f:
      names_from_f = [x.strip() for x in f]
   names,shell=loadXYZ(args.shell)
   
   shell_center = shell.mean(axis=0)
   shell-=shell_center
   
   np.random.shuffle(shell)
   
   names,mono=loadXYZ(args.monomer)
   
   n=args.n0
   
   mono = align_along_longest(mono,n)

   a=mono[n]
   
   residues = []
   
   V=args.scale*6.25/0.5291772
   dmin=args.d0
   s=[np.array([0.0,0.0,0.0])]
   minus=[0]
   N=0
   res, at = 0, 0
   with open(args.out,'w') as f:
      #f.write('%s\n\n'%(' '*72))
      for x in shell:
         rot,trans=get_rot(a,x,0)
         new=rot.dot(mono.T).T+trans+shell_center
         new_s=new[n]
         at=write(f, names, names_from_f, new, res, at, 'LIG')
         res+=1
         residues.append(new)
      #f.seek(0)
      #f.write('%i'%N)
   
   s= list(names).index('S')
   
   residues=np.array([[x[s]] for x in residues])
   dav, dmin = stats_of_two_arrs(residues)
   print 'Mono 1: ',dmin, dav
