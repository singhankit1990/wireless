classdef SISOMPChan < matlab.System
    % SISOMPChan:  SISO multi-path fading channel    
    properties 
        fsamp;   % Sample rate in Hz
                
        % Path properties
        gain;  % path gains in dB
        dly;   % delays in seconds
        dop;   % doppler shift of each path in Hz
        
        % Fractional delay object
        fracDly;
        
        % Initial set of phases for the next step call
        phaseInit;
                                
    end
    
    methods 
        function obj = SISOMPChan(varargin)
            % Constructor:  
         
            
            if nargin >= 1
                obj.set(varargin{:});
            end
            
        end
        
    end
    methods (Access = protected)
        function setupImpl(obj)
              
              obj.fracDly = dsp.VariableFractionalDelay("InterpolationMethod","Farrow","FilterLength",8,"FarrowSmallDelayAction","Use off-centered kernel","MaximumDelay",1024);

        end
        
        function resetImpl(obj)
           
            reset(obj.fracDly)
            
           %  Initializing  phases, phaseInit, to a row vector of 
            % dimension equal to the number of paths with uniform values 
            % from 0 to 2pi
            npath = length(obj.dly);
            obj.phaseInit = 2*pi*randn(npath,1);
        end
        
        function releaseImpl(obj)
            
           % Releasing  the fracDly object
            release(obj.fracDly);
        end
        
        function y = stepImpl(obj, x)
            % step:  Run a vector of samples through the channel
            
                        
            %  Computing  the delay in samples
                 dlySamp = obj.fsamp*obj.dly;
            
            %  Computing gain of each path in linear scale
                gainLin = db2pow(obj.gain);
             
            % delayed version 
                 xdly = obj.fracDly(x,dlySamp);
      
            
            % computing phase rotations 
                nsamp = length(x);
                npath = 24;
                phase = zeros((nsamp+1),npath);
                
                
               for i = 1 : npath
                   for j = 1: nsamp+1
                        phase(j,i) = -2*pi*obj.dop(i)*obj.dly(i);
                   end
               end

            
            % Saving  the final phase, phase(nsamp+1,:)
            % as phaseInit for the next step.

            obj.phaseInit = phase(nsamp+1,:);
            
            % Applying  the phases and gain to each path, adding  the
            % results
            y = zeros(nsamp,1);
            
            for k = 1:nsamp
                m =0;
                for j = 1: npath
               
                m = m + sqrt(gainLin(1,j))*exp(phase(1,j))*xdly(k,j);
                end
                y(k,1) = m;
            end
             
        end
    end
end