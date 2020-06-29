import argparse
parser=argparse.ArgumentParser(description='Make SAM', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--shell',type=str,help='XYZ of outer atoms')
parser.add_argument('--monomers', nargs='+', type=str, help='XYZ of monomer')
parser.add_argument('--monomers_resnames', nargs='+', type=str )
parser.add_argument('--monomers_numbers', nargs='+', type=int )
parser.add_argument('--backups', nargs=2, type=str, default=[] )
parser.add_argument('--names', nargs='+', type=str, help='names of monomer')
parser.add_argument('--anchor_idxes', type=int, nargs='+', help="index of monomer's atom to be attached", default=0)
parser.add_argument('--use_center', action='store_true')
parser.add_argument('--out', type=str, help="output file", default='out.pdb')

args=parser.parse_args()


from thomson import  gen_attachement_points
from put_on_sphere import put_all_on_shell
from sys import stdout
from simple_shell import line_fmt, attach_monomers, loadXYZ
import numpy as np

#============= READ DATA ================
print 'Shell config temorarly fixed'
shell_config = dict(resname='COR', atname='AU', sep=3.5, name='Au')
monomer_indices = args.anchor_idxes
monomer_numbers = args.monomers_numbers
number_of_points = sum(monomer_numbers)
Nmonomers = len(monomer_numbers)
assert Nmonomers<=2, 'so far up to two monomers supported'
ratio = float(monomer_numbers[0])/number_of_points
radius = 22
offset=0

_, shell = loadXYZ(args.shell)
monomer_resnames = args.monomers_resnames
monomer_names, monomer_elements, monomers = [], [], []

for i in range(Nmonomers):
   names, x_arr = loadXYZ(args.monomers[i])
   monomers.append(x_arr)
   monomer_elements.append(names)
   with open(args.names[i], 'r') as f:
      monomer_names.append([x.strip() for x in f])
   
print 'data read'
#============== get points ==========================
if args.backups ==[]:
   initial_guess = np.random.normal(0,1,size=(number_of_points, 3))
   initial_guess = initial_guess*radius/np.linalg.norm(initial_guess, axis=1).reshape(-1,1)
   qs=np.ones(initial_guess.shape[0])
   
   attachement_config = dict(enhance_half=True, fix_half=True, pull=True, maxit=10000, radius=radius, factor=ratio)
   qs, attach_xyz = gen_attachement_points(qs, initial_guess, **attachement_config)
   
   #put points on the shell
   attach_xyz = put_all_on_shell(attach_xyz, shell, shell_config['sep'])
   
   np.save('qs_backup.npy', qs)
   np.save('attach_xyz_backup.npy', attach_xyz)
else:
   qs = np.load(args.backups[0])
   attach_xyz = np.load(args.backups[1])

print 'attachement done'
#============= WRITE FILE =========================
resid=1
atid=1

if Nmonomers==2:
   idents = np.unique(qs)
   indices = [np.where(qs==i)[0] for i in idents]
elif Nmonomers==1:
   indices = [np.arange(len(qs))]
else:
   raise ValueError('At least one monomer should be present')

with open(args.out, 'w') as f:
   #write center
   if args.use_center:
      x,y,z = shell.mean(axis=0)   
      f.write(line_fmt%(atid, 'AC', 'CNT', resid, x, y, z, 'Au'))
      resid+=1
      atid+=1

   #write_shell
   for x, y, z in shell:
      f.write(line_fmt%(atid, shell_config['atname'], shell_config['resname'], resid, x, y, z, shell_config['name']))
      resid+=1
      atid+=1
   print 'shell written'

   #attach monomers
   for i in range(Nmonomers):
      attachement_points = attach_xyz[indices[i]]
      residues = attach_monomers(monomers[i], monomer_indices[i], attachement_points, offset)
      for crd in residues:
         for idx, (x,y,z) in enumerate(crd):
            f.write(line_fmt%(atid, monomer_names[i][idx], monomer_resnames[i], resid, x, y, z, monomer_elements[i][idx]))            
            atid+=1
         resid+=1
         stdout.write('written %i residues\r'%resid)

   print '\nDone'
      
         



