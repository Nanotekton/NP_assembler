import numpy as np

def get_rot(start,end,offset=0):
	rs=np.sqrt(start.dot(start))
	re=np.sqrt(end.dot(end))
	cos=start.dot(end)/(rs*re)
	os=np.cross(start,end)
	sin=np.sqrt(os.dot(os))/(rs*re)
	if sin!=0:
		os/=sin*rs*re
	else: 
		os=np.random.random(3)
		rzut=os.dot(start)
		os=os-rzut*start/rs**2
		mod=np.sqrt(os.dot(os))
		os/=mod
	sin=np.sqrt((1-cos)/2)
	cos=np.sqrt((1+cos)/2)
	
# unit quaternion ; teta-> rotation angle ;ux,uy,uz -> unit vector of rotation axis (left-handed)
# i;j;k -> imaginary units
# q= cos teta/2 + sin teta/2 (ux*i + uj*j +uk*k)
# q= (w,  x,  y,  z)
#
#Conversion to rotation matrix
# WATCH THE MATRIX!!!!
#
# |  1-2*(y**2 + z**2)		2*(x*y - z*w)		2*(x*z + y*w)       |
# |  2*(x*y + z*w)		1-2*(x**2 + z**2)	2*(y*z - x*w)       | 
# |  2*(x*z - y*w)		2*(y*z + x*w)		1-2*(y**2 + x**2)   |
#
	w=cos
	x,y,z=sin*os	  
	rot=np.array([[1-2*(y**2 + z**2), 2*(x*y - z*w), 2*(x*z + y*w)],[2*(x*y + z*w), 1-2*(x**2 + z**2),  2*(y*z - x*w)],[2*(x*z - y*w), 2*(y*z + x*w), 1-2*(y**2 + x**2)]])
	return rot,(re-rs+offset)*end/re

if __name__=='__main__':
	x=np.array([1.0,0.0,0.0])
	y=np.array([-1.0,0.0,0.0])
	rot,trans=get_rot(x,y)
	xp=rot.dot(x)+trans
	print xp
