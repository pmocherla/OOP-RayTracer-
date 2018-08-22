"""
@author: Priyanka Mocherla
@version: Python 3
@created: 12/02/2017

A ray tracer module to model simple optical systems using OOP. Primarily deals with
refracting surfaces. Example use of this module can be found by running this script.

To use this class the following libararies are required:
    - Numpy
    - matplotlib

"""
import numpy as np
import matplotlib.pyplot as plt

class Ray:
    """ Class to model a ray of light using position and direction vectors to trace its path"""
    
    def __init__(self, point, direction):
        """Initial ray position and direction
        
        Args:
            point (array-like) - Start point of ray
            direction (array-like) - Initial direction of ray
        """
        self._p = np.array([point])
        self._k = np.array([direction])
        
    def p(self):
        """Returns the current position of the ray"""
        return self._p[-1]
        
    def k(self):
        """Returns the current direction of the ray"""
        return self._k[-1]
        
    def append(self, p, k):
        """Adds additional positions and directions to the ray
        
        Args:
            p (array-like) - Position of ray
            k (array-like) - Direction of ray
        """
        self._p = np.append(self._p, [p] , axis=0)
        self._k = np.append(self._k, [k], axis=0)
        
    def vertices(self):
        """Returns the full traced path taken by the ray"""
        return self._p
        

class OpticalElement:
    """Class to model optical elements such as lenses"""
    
    def propagate_ray(self, ray):        
        "propagate a ray through the optical element"
        raise NotImplementedError()
        
      
class SphericalRefraction(OpticalElement):
    """Class to model a spherical refracting object"""
    
    def __init__(self, z0, curvature, n1, n2, aperture_rad):
        """ Initialise lens parameters
        
        Args:
            z0 (float) - intercept of the surface with the z axis
            curvature (float) - curvature of the surface, equal to the inverse of the radius of curvature.
                                positive for convex edges, negative for concave.
            n1 (float) - refractive index left side of the surface
            n2 (float) - refractive index right side of the surface
            aperture_rad - vertical extent of lens
        """                      
        self._z0 = z0
        self._curvature = curvature
        self._n1 = n1
        self._n2 = n2
        self._aperture_rad = aperture_rad
        
    def _intercept(self, ray):
        """ Returns intercept of the ray with the lens, if any.
        
        Args:
            ray (instance) - instance of the Ray class
            
        Returns:
            (array-like) position vector of the intersection of the ray and the lens, if any.
            
        Example:
            a = Ray([0,0,0], [0,0,1])
            b = SphericalRefraction(1,1,1,1,1)
            print(b.intercept(a))
        """
        #Calculate the normal of the incident direction vector
        self._k_hat = ray._k[-1] / np.linalg.norm(ray._k[-1])
        
        #Determine whether the lens edge is convex (1) or concave (-1) from left to right
        if self._curvature < 0:
            const = -1
        else:
            const = 1
        
        #Check for curved surface - plane surfaces are not applicable in the following calculations
        if self._curvature != 0:
            
            #Calculate radius of lens and distance from inital ray point to centre of lens
            self._radius = 1/float(self._curvature)
            self._r = ray._p[-1] - np.array([0, 0, self._radius]) - np.array([0 ,0, self._z0])
            dot = np.dot(self._r, self._k_hat)
            sqr = dot*dot - (np.dot(self._r, self._r) - self._radius*self._radius)
            
            #Check whether a valid intercept is possible
            if sqr >= 0 :
                length = -1*dot - const * np.sqrt(sqr)
            
                #Calculate intercept
                self._inter = ray._p[-1] + length*self._k_hat
                
                #Determine whether aperture rad is greater than the calculated intercept
                if abs(self._inter[0]) <= self._aperture_rad:
                    return self._inter
                
                else:
                    return None
                
            else: 
                return None
            
        else:
            #Calculate intercept for plane surfaces
            length = self._z0 - ray._p[-1][-1]
            self._inter = ray._p[-1] + length/float(ray._k[-1][-1]) * ray._k[-1]
            
            #Determine whether aperture rad is greater than the calculated intercept
            if abs(self._inter[0]) < self._aperture_rad:
                return self._inter
            else:
                return None

            
    def propagate_ray(self, ray):
        """Propagate ray through spherical object
        
        Args:
            ray (instance) - instance of the Ray class
            
        Example:
            a = Ray([0,0,0], [0,0,1])
            b = SphericalRefraction(1,1,1,1,1)
            print(b.propagate_ray(a))
        """
        #Calculate the intercept and check if it is valid
        intercept = self._intercept(ray)
        if intercept is not None:
            #Determine plane or curved surface normal vector
            if self._curvature != 0:
                normal = (intercept - (np.array([0,0,self._z0]) + np.array([0, 0, self._radius])))
            else:
                normal = np.array([0,0,-1])
            
            #Calculate output direction of ray after propagating through lens
            out_k = snell(ray._k[-1], normal, self._n1, self._n2)
            
            #If valid output vector, append position and direction to ray instance
            if out_k is not None:
                ray.append(intercept, out_k)
                
                    
                        
class OutputPlane(OpticalElement):
        """Class modelling an output plane to measure intersection of rays with plane.
           For visualisation purposes."""
        
        def __init__(self, z):
            """ Initialise output plane parameters
        
            Args:
                z (float) - intercept of the plane with the z axis
            """
            self._z = z
            
        def _intercept(self, ray):
            """ Returns intercept of the ray with the output plane.
        
            Args:
                ray (instance) - instance of the Ray class
            
            Returns:
                (array-like) position vector of the intersection of the ray and the output plane.
            
            Example:
                a = Ray([0,0,0], [0,0,1])
                b = OutputPlane(19)
                print(b.intercept(a))
            """
            length = self._z - ray._p[-1][-1]
            self._inter = ray._p[-1] + length/float(ray._k[-1][-1]) * ray._k[-1]
            return self._inter
        
        def propagate_ray(self, ray):
            """ Propagates ray to output plane.
        
            Args:
                ray (instance) - instance of the Ray class
            
            Example:
                a = Ray([0,0,0], [0,0,1])
                b = OutputPlane(19)
                print(b.propagate_ray(a))
            """
            ray.append(self._intercept(ray), ray._k[-1])
            
            
            

def snell(incident_k, surface_normal, n1, n2):
    """3D Snells law calculation of output direction of incident vector through a surface.
    
    Args:
        inident_k (array-like) - incident ray direction vector
        surface_normal (array-like) - vector normal to refracting surface
        n1 (float) - refractive index on left side of surface
        n2 (float) - refractive index on right side of surface
        
    Returns:
        (array-like) output direction vector through the surface
    """
    #Calculate the angle between the ray and the surface normal
    theta = np.arccos(np.dot(incident_k, surface_normal) / (np.linalg.norm(surface_normal) * np.linalg.norm(incident_k)))
    n2 = float(n2)
    
    #Condition for total internal reflection must not be satisfied.
    if np.sin(theta) < n2/n1:
        #Normalise the vectors
        inc_k_hat = incident_k / np.linalg.norm(incident_k) 
        surf_norm_hat =  surface_normal / np.linalg.norm(surface_normal)

        #Calculation of Snells law broken into two parts and combined
        one = n1/n2* np.cross(surf_norm_hat, np.cross(surf_norm_hat, inc_k_hat))
        two = surf_norm_hat * np.sqrt(1 - (n1/n2)**2 * (np.dot(np.cross(surf_norm_hat, inc_k_hat),np.cross(surf_norm_hat, inc_k_hat))))
        return one - two
        
    else:
        return None
        
def rtpairs(R, N):
    """Generator to produce polar coordinates
    
    Args:
        R (list) - radii from origin for polar coordinates
        N (list) - number of points to plot at each radius
        
    Yields: 
        tuple: radius of coordinate and list of phases at specified radius
 
    Examples: 
        R = [0.0, 0.1, 0.2]
        T = [1, 10, 20]
        for r, t in genpolar.rtpairs(R, T):
            plt.plot(r * cos(t), r * sin(t), 'bo')   
    """
    for i in range(len(R)):
        r = R[i]
        t = (360./N[i] * np.pi/180.)*np.arange(1,N[i]+ 1)
        yield r, t

def rtuniform(n, rmax, m):
    """Produces uniform polar coordinates
    
    Args:
        n (int) - number of rings to plot between maximum radius and origin
        rmax (float) - maximum radius of polar coordinates
        m (int) - density of points at each radius increases by m each step from origin
        
    Returns: 
        tuple: radius of coordinate and list of phases at specified radius
        
    Examples:
        for r, t in rtuniform(10, 10, 6):
        plt.plot(r * np.cos(t), r * np.sin(t), 'bo')
    """
    ri = rmax/float(n) *np.arange(0,n+1)
    ni = np.concatenate([[1], m* np.arange(1, n+1)])
    return rtpairs(ri,ni)


#------------------------- Example Usage ------------------------#      
def main():
    #Modelling the diffracion of a plano convex lens, curved surface facing right
    b = SphericalRefraction(100, -0.03, 1.5, 1.0, 30) #Lens placed at 100mm
    c = OutputPlane(163) #Placed at focal point
    f, (ax1, ax2) = plt.subplots(1, 2)
    for r, t in rtuniform(5, 10, 5):
        for i in range(len(t)):
            a = Ray([r * np.cos(t[i]), r * np.sin(t[i]), 0], [0, 0.00, 1])
            b.propagate_ray(a)
            c.propagate_ray(a)
            
            ax1.plot(a.vertices()[:,2], a.vertices()[:,0], 'r')
            ax2.plot(a.vertices()[-1][1], a.vertices()[-1][0], 'ro')
            
    ax1.set_title("Ray path")
    ax1.set_ylabel("x / mm")
    ax1.set_xlabel("z / mm")
    ax2.set_title("Output plane (Focus)")
    ax2.set_ylabel("x / mm")
    ax2.set_xlabel("y / mm")
    
    plt.show()


if __name__ == "__main__":
    main()

