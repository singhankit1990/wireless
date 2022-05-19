classdef ElemWithAxes < matlab.System
    % ElemWithAxes:  An antenna element with a local frame of reference
   
    properties
        % The antenna object from the phased array toolbox
        ant = [];
        
        % Azimuth and elevation angle of the element peak directivity
        axesAz = 0;
        axesEl = 0;
        
        % Axes of the element local coordinate frame of reference
        axesLoc = eye(3);
        
        % Frequency in Hz
        fc = 0;
        
        % Directivity interpolant
        dirInterp = [];
        
        % Velocity vector in 3D in m/s
        vel = zeros(1,3);
        
    end
    
    methods
        function myobj = ElemWithAxes(fc, ant)
            % Constructor 
          
                myobj.ant = ant;
                myobj.fc = fc;
        end
        
        function alignAxes(myobj,az,el)
            % Aligns the axes to given az and el angles
            
         
            myobj.axesAz = az;
            myobj.axesEl = el;
            
            
          
            myobj.axesLoc = azelaxes(az,el);
            
        end
        
        function dop = doppler(myobj,az,el)
            % Computes the Doppler shift of a set of paths 
            % The angles of the paths are given as (az,el) pairs
            % in the global frame of reference.
            
         
     
            [x,y,z] = sph2cart(deg2rad(az),deg2rad(el),1);
            u = [x;y;z];
        
            myobj.vel  = myobj.vel;
            vp = physconst('lightspeed');
            for i = 1 : length(u)
                dop(i) = myobj.vel*u(:,i).*myobj.fc/vp;
            % where vc = speed of light
            
            
            end
        
        end
    end
    
    methods (Access = protected)
         function setupImpl(myobj)
            % setup:  This is called before the first step.
            
            
            % Getting the pattern from ant.pattern
            [direct,az,el] = myobj.ant.pattern(myobj.fc,"Type","directivity");
            
            % Creating the gridded interpolant object. 
            
            
            myobj.dirInterp = griddedInterpolant({el,az},direct);
            
            
            
            
        end
        function dir = stepImpl(myobj, az, el)
            
           
          
            
            .
            ug = [az;el;ones(1,length(az))];
            Loc = global2localcoord(ug,'ss',zeros(3,1),myobj.axesLoc);
            azLoc = Loc(1,:);
            elLoc = Loc(1,:);
            
       
            dir = myobj.dirInterp(azLoc,elLoc);
        end
        
    end
    
end

