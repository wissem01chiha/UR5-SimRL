import numpy as np
import transform

class Robot:
    """
    Container class of the robot parametric description.
    """
    def __init__(self, name, nl, nj , directory= None) -> None:

        # name of the robot: string 
        self.name = name
        # number of links """
        self.nl = nl
        """ directory name"""
        self.directory = directory
        # gravity vector: 3x1 matrix 
        self.G = np.array([0, 0, -9.81])
        self.nj = nj
        # links dh paramters 
        self.a = list(None for j in range(nj))
        self.d = list(None for j in range(nj))
        self.alpha = list(None for j in range(nj))
        # link inertia tensors
        self.inertia_mats = list(None for j in range(nj))
        # link mass values 
        self.mass = list(None for j in range(nj))
        # joint positions
        self.q = list(None for j in range(nj))
        # joint velocities 
        self.qdot = list(None for j in range(nj))
        # joint accelerations
        self.qddot = list(None for j in range(nj))
        # link velocities at com
        self.vels = list(None for j in range(nj))
        # link accelerations at com 
        self.accels = list(None for j in range(nj))
        # joint corlois torques 
        self.corils = list(None for j in range(nj))
        # joint gravitational torques
        self.grv_torques = list(None for j in range(nj))
        # joint applied torques
        self.torque = list(None for j in range(nj))
        # joint friction torques
        self.frt_torque = list(None for j in range(nj))
        # joint stifness coefficents
        self.stiffness = list(None for j in range(nj))
        
        # mass matrix 
        self.mass_mat = np.zeros((nj,nj))
        # corlois matrix
        self.corils_mat = np.zeros((nj,nj))
        
        
    def __str__(self):
        str_format = ""
        str_format = str_format + "Robot ({0}):\n".format( 
            self.model_type )
            
        str_format = str_format + "-------------------\n"
        attrs = [
            attr for attr in dir(self) \
            if not attr.startswith('_')
        ]
        for attr in attrs:
            items = getattr(self, attr)
            if hasattr(items, '__iter__'):
                attr_str = self._str_items(items)
            else:
                attr_str = '\t' + str(items) + '\n'
            str_format = str_format + str(attr) + ": \n"
            str_format = str_format + attr_str
        return str_format
    
    def _str_items(self, items):
            """
            Create a string representation for a given attribute.

            Args:
                items: An attribute of the class which is a list

            Returns:
                A string representation of the attribute.
            """
            row_format = '\t' + ('{0:^6} : {1}') + '\n'
            str_format = ""
            for idx, item in enumerate(items):
                str_format = str_format + row_format.format(*(
                    str(idx), str(item)
                ))
            return str_format  
        
    @classmethod      
    def rotation_matrix(self, i):
        """
        
        """
        return 
    
    @classmethod
    def homogenus_matrix(self, i):
        """
        
        """
        return
    
    @classmethod
    def translation_vector(self,i):
        
        return 
    
    @classmethod
    def link_com_vect(self,i):
        """
        Returns the link "i" center of mass coordiantes wrt
        joint articulated frame. 
        """
        return 
    
    @classmethod
    def link_inertia_mat(self,i ):
        """ """
        return 
    